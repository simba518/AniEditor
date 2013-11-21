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
void function1_grad(const real_1d_array &x, double &func, real_1d_array &grad,void*ptr) {

  func = 100*pow(x[0]+3,4) + pow(x[1]-3,4);
  grad[0] = 400*pow(x[0]+3,3);
  grad[1] = 4*pow(x[1]-3,3);
}


bool VolApproxTri::generateKeyVolMesh(){
  
  real_1d_array x = "[0,0]";
  double epsg = 0.0000000001;
  double epsf = 0;
  double epsx = 0;
  ae_int_t maxits = 0;
  minlbfgsstate state;
  minlbfgsreport rep;

  minlbfgscreate(1, x, state);
  minlbfgssetcond(state, epsg, epsf, epsx, maxits);
  alglib::minlbfgsoptimize(state, function1_grad);
  minlbfgsresults(state, x, rep);

  printf("%d\n", int(rep.terminationtype)); // EXPECTED: 4
  printf("%s\n", x.tostring(2).c_str()); // EXPECTED: [-3,3]

  return false;
}
