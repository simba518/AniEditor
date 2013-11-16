#include <vector>
#include <MatrixTools.h>
#include <assertext.h>
#include <AuxTools.h>
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

  Matrix3d F,Q,S;
  createFromRowMajor(F,m);
  ModifiedPD3x3 (F,Q,S);

  // symmetric part
  y[3+0] = S(0,0)-1.0f;
  y[3+1] = S(0,1);
  y[3+2] = S(0,2);
  y[3+3] = S(1,1)-1.0f;
  y[3+4] = S(1,2);
  y[3+5] = S(2,2)-1.0f;

  // antisymmetric part
  // http://en.wikipedia.org/wiki/Rotation_matrix#Quaternion
#define COPYSIGN(x,y) (y>0?1:-1)*fabs(x)

  const double t = Q(0,0)+Q(1,1)+Q(2,2); assert_ge(1.0f+t,0.0f);
  const double r = sqrt(1.0f+t);
  double qw= 0.5*r;
  qw = (qw > 1.0) ? 1.0:qw; /// @todo need more elegant method to compute the Quaternion.
  qw = (qw < -1.0) ? -1.0:qw;
  const double qx= COPYSIGN(0.5*sqrt(fabs(1+Q(0,0)-Q(1,1)-Q(2,2))),Q(2,1)-Q(1,2));
  const double qy= COPYSIGN(0.5*sqrt(fabs(1-Q(0,0)+Q(1,1)-Q(2,2))),Q(0,2)-Q(2,0));
  const double qz= COPYSIGN(0.5*sqrt(fabs(1-Q(0,0)-Q(1,1)+Q(2,2))),Q(1,0)-Q(0,1));
  assert_in((qw*qw+qx*qx+qy*qy+qz*qz),1.0f-1e-10,1.0f+1e-10);

  // http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/
  const double temp = sqrt(1.0-qw*qw);
  const double theta = 2.0f*acos(qw);
  assert_eq(theta,theta);
  if (temp > 1e-10){ /// @todo need more elegant method to compute the Quaternion.
  	y[0] = theta*qx/temp;
  	y[1] = theta*qy/temp;
  	y[2] = theta*qz/temp;
  }else{
	y[0] = theta*qx;
  	y[1] = theta*qy;
  	y[2] = theta*qz;
  }
}

void RSCoordComp::RS2NodeRotVec(pTetMesh_const tetmesh,const VectorXd &y, VectorV3 &w){
  
  assert(tetmesh != NULL);
  assert_eq(tetmesh->tets().size(),(int)y.size()/9);

  w.resize(tetmesh->nodes().size());
  vector<int> eleCount(tetmesh->nodes().size());
  for (size_t i = 0; i < w.size(); ++i){
	w[i].setZero();
	eleCount[i] = 0;
  }
	  
  for (int e = 0; e < tetmesh->tets().size(); ++e){
	for (int v = 0; v < 4; ++v){
	  const int i = tetmesh->tets()[e][v];
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

void RSCoordComp::RS2NodeRSVec(pTetMesh_const tetmesh,const VectorXd &y, VectorV6 &s){
  
  assert(tetmesh != NULL);
  assert_eq(tetmesh->tets().size(),(int)y.size()/9);

  s.resize(tetmesh->nodes().size());
  vector<int> eleCount(tetmesh->nodes().size());
  for (size_t i = 0; i < s.size(); ++i){
	s[i].setZero();
	eleCount[i] = 0;
  }
	  
  for (int e = 0; e < tetmesh->tets().size(); ++e){
	for (int v = 0; v < 4; ++v){
	  const int i = tetmesh->tets()[e][v];
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
