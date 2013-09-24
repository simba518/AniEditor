#include <vector>
#include <MatrixTools.h>
#include <assertext.h>
#include <AuxTools.h>
#include <polarDecomposition.h>
#include "RSCoordComp.h"
using namespace std;
using namespace UTILITY;
using namespace LSW_WARPING;
using namespace EIGEN3EXT;

void RSCoordComp::constructWithWarp(const SparseMatrix<double> &G, const VectorXd &p, VectorXd &y){
  
  assert_eq(G.cols(), p.size());
  assert_eq(G.rows()%9, 0);

  const VectorXd m = G*p;
  y.resize(m.size());
  const int elem_num = m.size()/9;
  for (int j =0; j < elem_num; ++j){
    const double *mj = &(m[9*j]);
	double *yj = &(y[9*j]);
	constructWithWarp_OneTet(mj, yj);
  }
}

void RSCoordComp::constructWithoutWarp(const SparseMatrix<double> &G, const VectorXd &u, VectorXd &y){
  
  const int elem_num = G.rows()/9;
  assert_gt(elem_num,0);
  assert_eq(G.rows()%9, 0);
  assert_eq(G.cols(),u.size());

  VectorXd m = G*u;
  for (int i = 0; i < elem_num; ++i){
    m[i*9 + 0] += 1.0;
    m[i*9 + 4] += 1.0;
    m[i*9 + 8] += 1.0;
  }
  y.resize(elem_num*9);
  for (int j =0; j < elem_num; ++j){
    const double *mj = &(m[9*j]);
	double *yj = &(y[9*j]);
	constructWithoutWarp_OneTet(mj, yj);
  }
}

void RSCoordComp::constructWithWarp_OneTet(const double *m, double *y){

  // rotation part
  const double w0 = (m[7]-m[5])/2;
  const double w1 = (m[2]-m[6])/2;
  const double w2 = (m[3]-m[1])/2;
  y[0] = w0;   y[1] = w1;   y[2] = w2;
  
  // strain part
  const double s0 = m[0];
  const double s3 = m[4];
  const double s5 = m[8];
  const double s1 = (m[1]+m[3])/2;
  const double s2 = (m[2]+m[6])/2;
  const double s4 = (m[5]+m[7])/2;
  y[3] = s0;   y[4] = s1;   y[5] = s2;
  y[6] = s3;   y[7] = s4;   y[8] = s5;

}

void RSCoordComp::constructWithoutWarp_OneTet(const double *m, double *y){

  Matrix3d Fm,Qm,Sm;
  createFromRowMajor(Fm,m);
  ModifiedPD3x3 (Fm,Qm,Sm);

  double Q[9], S[9]; // row major 3x3
  createRowMajor(Qm,Q);
  createRowMajor(Sm,S);

  assert_in (Qm.determinant(),1.0f-1e-3,1.0f + 1e-3);
  assert_lt ((Qm.transpose()*Qm - Matrix3d::Identity()).norm(), 1e-3);
  assert_lt ((Sm.transpose()-Sm).norm(), 1e-3);
  assert_lt ( (Fm-Qm*Sm).norm(), 1e-2 );

  // antisymmetric part
  const double qw = sqrt(1 + Q[0] + Q[4] + Q[8]) / 2.0;
  const double qx = (Q[7] - Q[5])/( 4 *qw);
  const double qy = (Q[2] - Q[6])/( 4 *qw);
  const double qz = (Q[3] - Q[1])/( 4 *qw);
  
  const double theta = 2*acos(qw);
  if (fabs(theta) > 1e-13){

	const double temp = sqrt(1.0 - qw*qw);
	assert_gt(temp,1e-13);
	const double _x = qx/temp;
	const double _y = qy/temp;
	const double _z = qz/temp;

	const double _norm = sqrt(_x*_x + _y*_y + _z*_z);
	assert_gt(_norm,1e-13);
	y[0] = theta * _x /_norm;
	y[1] = theta * _y /_norm;
	y[2] = theta * _z /_norm;
  }else{
	y[0] = 0.0;
	y[1] = 0.0;
	y[2] = 0.0;	
  }
  
  // symmetric part
  y [3+0] = S[0] - 1.0f;
  y [3+1] = S[1];
  y [3+2] = S[2];
  y [3+3] = S[4] - 1.0f;
  y [3+4] = S[5];
  y [3+5] = S[8] - 1.0f;
}

void RSCoordComp::RS2NodeRotVec(pVolumetricMesh_const tetmesh,const VectorXd &y, VectorV3 &w){
  
  assert(tetmesh != NULL);
  assert_eq(tetmesh->numElements(),(int)y.size()/9);

  w.resize(tetmesh->numVertices());
  vector<int> eleCount(tetmesh->numVertices());
  for (size_t i = 0; i < w.size(); ++i){
	w[i].setZero();
	eleCount[i] = 0;
  }
	  
  for (int e = 0; e < tetmesh->numElements(); ++e){
	for (int v = 0; v < 4; ++v){
	  const int i = tetmesh->vertexIndex(e,v);
	  assert_in(i,0,(int)w.size()-1);
	  w[i] += y.segment(e*9,3);
	  eleCount[i] ++;
	}
  }

  for (size_t i = 0; i < w.size(); ++i){
	assert_gt (eleCount[i],0);
	w[i] /= (double)eleCount[i];
  }
}

void RSCoordComp::RS2NodeRSVec(pVolumetricMesh_const tetmesh,const VectorXd &y, VectorV6 &s){
  
  assert(tetmesh != NULL);
  assert_eq(tetmesh->numElements(),(int)y.size()/9);

  s.resize(tetmesh->numVertices());
  vector<int> eleCount(tetmesh->numVertices());
  for (size_t i = 0; i < s.size(); ++i){
	s[i].setZero();
	eleCount[i] = 0;
  }
	  
  for (int e = 0; e < tetmesh->numElements(); ++e){
	for (int v = 0; v < 4; ++v){
	  const int i = tetmesh->vertexIndex(e,v);
	  assert_in(i,0,(int)s.size()-1);
	  s[i] += y.segment(e*9+3,6);
	  eleCount[i] ++;
	}
  }

  for (size_t i = 0; i < s.size(); ++i){
	assert_gt (eleCount[i],0);
	s[i] /= (double)eleCount[i];
  }
}
