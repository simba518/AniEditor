#include <alglib/stdafx.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <alglib/optimization.h>
#include "VolApproxTri.h"
#include <MatrixIO.h>
using namespace UTILITY;
using namespace alglib;

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

bool VolApproxTri::generateKeyVolMesh(pObjmesh_const objMeshKey){

  // precompute
  this->prepare(objMeshKey);

  pTetMesh_const tetMeshRest = _embed->getTetMesh();
  pObjmesh_const objMeshRest = _embed->getObjMesh();

  // set optimization parameters.
  double epsg = 1e-12;
  double epsf = 0.0f;
  double epsx = 0.0f;
  ae_int_t maxits = 0;
  minlbfgsstate state;
  minlbfgsreport rep;

  const int n = tetMeshRest->nodes().size();
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

  return alglibErrorReport(int(rep.terminationtype));
}

void VolApproxTri::prepare(pObjmesh_const objMeshKey){

  pTetMesh_const tetMeshRest = _embed->getTetMesh();
  pObjmesh_const objMeshRest = _embed->getObjMesh();
  _elasticModel = pElasticForceTetFullStVK(new ElasticForceTetFullStVK(tetMeshRest));
  tetMeshRest->buildInterpMatrix(_embed->getInterpNodes(),_embed->getInterpWeights(),tetMeshRest->nodes().size(),_A);
  _objKeyU = objMeshKey->getVerts()-objMeshRest->getVerts();
  tetMeshRest->nodes(_volRestU);
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

bool VolApproxTri::alglibErrorReport(const int code)const{
  
  switch(code){
  case -7:cout << "gradient verification failed.\n"; break;
  case -3:cout << "inconsistent constraints.\n"; break;
  case 1:cout << "relative function improvement is no more than EpsF.\n"; break;
  case 2:cout << "relative step is no more than EpsX.\n"; break;
  case 4:cout << "gradient norm is no more than EpsG.\n"; break;
  case 5:cout << "MaxIts steps was taken.\n"; break;
  case 7:cout << "stopping conditions are too stringent,further improvement is impossible.\n"; break;
  default: cout << "the stop code is unkown: " << code << endl;
  }
  return code > 0;
}
