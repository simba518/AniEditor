#include <alglib/stdafx.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <alglib/optimization.h>
#include "VolApproxTri.h"
#include <MatrixIO.h>
using namespace UTILITY;
using namespace alglib;

bool VolApproxTri::loadRestVolMesh(const string filename){
  const bool succ = _tetMeshRest->load(filename);
  ERROR_LOG_COND("failed to load tetrahedron mesh from: "<<filename,succ);
  if (succ){
	_embed.setTetMesh(_tetMeshRest);
	_volKeyU.resize(_tetMeshRest->nodes().size()*3);
	_volKeyU.setZero();
  }
  return succ;
}

bool VolApproxTri::loadRestObjMesh(const string filename){
  const bool succ = _objMeshRest->load(filename);
  ERROR_LOG_COND("failed to load obj mesh from: "<<filename,succ);
  if (succ)
	_embed.setObjmesh(_objMeshRest);
  return succ;
}

bool VolApproxTri::loadKeyObjMesh(const string filename){
  const bool succ = _objMeshKey->load(filename);
  ERROR_LOG_COND("failed to load obj mesh from: "<<filename,succ);
  return succ;
}

bool VolApproxTri::loadInterpWeights(const string filename){
  
  const bool succ = _embed.loadWeights(filename);
  ERROR_LOG_COND("failed to load interpolation weights from:"<<filename,succ);
  return succ;
}

bool VolApproxTri::saveAll(const string filename){
  
  const bool succ1 = EIGEN3EXT::write(filename+".b",_volKeyU);
  ERROR_LOG_COND("failed to save the key u to: "<<filename,succ1);
  const bool succ2 = _tetMeshRest->writeVTK(filename+".vtk",_volKeyU);
  ERROR_LOG_COND("failed to save the key u to: "<<filename,succ2);
  return succ2&succ1;
}

// optimization using alglib
void function_grad(const real_1d_array &x, double &func, real_1d_array &grad,void*ptr) {

  VolApproxTri *appr = (VolApproxTri *)ptr;
  assert_gt(x.length(),0);
  const VectorXd &u = Map<VectorXd>((double*)&x[0],x.length());
  VectorXd g;
  appr->funGrad(u,func,g);
  assert_eq(grad.length(),g.size());
  memcpy(&grad[0],&g[0],sizeof(double)*g.size());
}

void function_func(const real_1d_array &x, double &func, void *ptr) {
  
  VolApproxTri *appr = (VolApproxTri *)ptr;
  assert_gt(x.length(),0);
  const VectorXd &u = Map<VectorXd>((double*)&x[0],x.length());
  appr->fun(u,func);
}

bool VolApproxTri::generateKeyVolMesh(){

  // precompute
  this->prepare();

  // set optimization parameters.
  double epsg = 1e-12;
  double epsf = 0.0f;
  double epsx = 0.0f;
  ae_int_t maxits = 0;
  minlbfgsstate state;
  minlbfgsreport rep;

  const int n = _tetMeshRest->nodes().size();
  real_1d_array u;
  u.setlength(n*3);
  for (int i = 0; i < n*3; ++i){
    u[i] = 0.0f;
  }

  // solve
  minlbfgscreate(1, u, state);
  minlbfgssetcond(state, epsg, epsf, epsx, maxits);
  alglib::minlbfgsoptimize(state, function_grad,NULL,(void*)this);
  minlbfgsresults(state, u, rep);

  // get results
  _volKeyU.resize(n*3);
  for (int i = 0; i < n*3; ++i){
    _volKeyU[i] = u[i];
  }

  printf("%d\n", int(rep.terminationtype)); 
  // printf("%s\n", u.tostring(2).c_str());

  return true;
}

bool VolApproxTri::generateKeyVolMeshNumDiff(){

  cout << "warnning: numerical differential method will be used, which is very slow!" << endl;
  
  // precompute
  this->prepare();

  // set optimization parameters.
  double epsg = 0.0000000001;
  double epsf = 0;
  double epsx = 0;
  double diffstep = 1.0e-6;
  ae_int_t maxits = 0;
  minlbfgsstate state;
  minlbfgsreport rep;

  const int n = _tetMeshRest->nodes().size();
  real_1d_array u;
  u.setlength(n*3);
  for (int i = 0; i < n*3; ++i){
    u[i] = 0.0f;
  }

  // solve
  minlbfgscreatef(1, u, diffstep, state);
  minlbfgssetcond(state, epsg, epsf, epsx, maxits);
  alglib::minlbfgsoptimize(state, function_func,NULL,(void*)this);
  minlbfgsresults(state, u, rep);

  // get results
  _volKeyU.resize(n*3);
  for (int i = 0; i < n*3; ++i){
    _volKeyU[i] = u[i];
  }

  printf("%d\n", int(rep.terminationtype)); 
  // printf("%s\n", u.tostring(2).c_str());
}

void VolApproxTri::funGrad(const VectorXd &u,double &fun,VectorXd &grad){

  // function
  const VectorXd du = (_A*u-_objKeyU);
  fun = 0.5f*du.dot(du);
  const VectorXd x = _volRestU+u;
  fun += _penalty*_elasticModel->energy(x);

  // grad
  _elasticModel->force(x,grad);
  grad = _A.transpose()*du - _penalty*grad;
}

void VolApproxTri::fun(const VectorXd &u,double &fun){

  // function
  const VectorXd du = (_A*u-_objKeyU);
  fun = 0.5f*du.dot(du);
  const VectorXd x = _volRestU+u;
  fun += _penalty*_elasticModel->energy(x);
}

void VolApproxTri::prepare(){
 
  _elasticModel = pElasticForceTetFullStVK(new ElasticForceTetFullStVK(_tetMeshRest));
  _tetMeshRest->buildInterpMatrix(_embed.getInterpNodes(),_embed.getInterpWeights(),_tetMeshRest->nodes().size(),_A);
  _objKeyU = _objMeshKey->getVerts()-_objMeshRest->getVerts();
  _tetMeshRest->nodes(_volRestU);
}
