#include <eigen3/Eigen/Dense>
#include <assertext.h>
#include <MatrixIO.h>
#include <ConMatrixTools.h>
#include "RS2Euler.h"
#include "DefGradOperator.h"
#include "ComputeBj.h"
using namespace Eigen;
using namespace LSW_WARPING;

/**
 * 1. compute the deformation gradient operator G
 * 2. compute the volumes V and Sqrt_V of all the tetrahedron.
 */
void RS2Euler::setTetMesh(pVolumetricMesh_const tetmesh){
  
  assert (tetmesh != NULL);
  this->tetmesh = tetmesh;
  const int elem_num = tetmesh->numElements();
  const int node_num = tetmesh->numVertices();
  SparseMatrix<double> V;

  if (DefGradOperator::compute(tetmesh,G)){

	assert_eq(G.cols(), node_num*3);
	assert_eq(G.rows(), elem_num*9);
	constructVolumeMatrix(tetmesh,V,Sqrt_V);
	assert_eq(V.cols(), G.rows());

	VG_t = (V*G).transpose();
	VG_t.makeCompressed();
	L = VG_t*VG_t.transpose();
  }else{
	ERROR_LOG("failed to compute gradient operator G");
  }
}

// compute the constraint matrix for the fixed nodes
void RS2Euler::setFixedNodes(const vector<int> &fixed_nodes){

  assert (tetmesh != NULL);
  const int node_num = tetmesh->numVertices();
  UTILITY::computeConM(fixed_nodes,F, node_num);
  assert_eq(F.cols(), node_num*3);
}

// compute the constraint matrix for the constrained nodes
void RS2Euler::setConNodes(const vector<set<int> > &c_nodes,const VectorXd &bcen_uc){

  assert (tetmesh != NULL);
  assert_eq ((int)c_nodes.size()*3, bcen_uc.size());
  this->barycenter_uc = bcen_uc;
  const int node_num = tetmesh->numVertices();
  if (c_nodes.size() > 0){
	UTILITY::computeBaryCenterConM(c_nodes,node_num,C);
  }else{
	C.resize(0,0);
  }
}

/**
 * assemble A and make factorization.
 */
bool RS2Euler::precompute(){
  
  bool succ = true;
  assembleA();
  A.makeCompressed();
  A_solver.compute(A);
  if(A_solver.info()!=Eigen::Success) {
	ERROR_LOG("failed to decompse A for solving A X = B");
	succ = false;
  }
  return succ;
}

/** 
 * reconstruct the displacements u in Euler space with RS coordinates y
 * provided.
 * 
 * @param y the RS coordinates.
 * @param u the constructed Euler coordinates represents the displacements.
 * 
 * @return true if construction is success.
 */
bool RS2Euler::reconstruct(const VectorXd &y, VectorXd &u){

  assert (tetmesh != NULL);
  const int node_numx3 = tetmesh->numVertices()*3;
  bool succ = true;

  // assemble the right_side
  VectorXd b;
  assemble_b(y,b);
  assert_eq(VG_t.cols(),b.size());
  assert_eq(VG_t.rows(),node_numx3);

  VectorXd right_side(node_numx3 + numFixedDofs());
  right_side.segment(0, node_numx3) = VG_t*b;
  right_side.segment(node_numx3,numFixedDofs()).setZero();
  right_side.segment(node_numx3,barycenter_uc.size()) = barycenter_uc;

  // solve A*u = right_side
  assert_eq (right_side.size(), A_solver.rows());
  u = A_solver.solve(right_side);

  // get the result 
  if(A_solver.info()!=Eigen::Success) {
	
  	succ = false;
  	u.resize(node_numx3);
  	u.setZero();
  	ERROR_LOG("failed to solve for A X = P.");
  }else{
	assert_gt(u.size(), node_numx3);
	const VectorXd x = u.head(node_numx3);
	u = x;
  }
  return succ;
}

void RS2Euler::constructVolumeMatrix(pVolumetricMesh_const tetmesh, SparseMatrix<double> &V,vector<double> &Sqrt_V)const{
  
  const int elem_num = tetmesh->numElements();

  Sqrt_V.resize(elem_num);
  V.resize(elem_num*9,elem_num*9);
  V.reserve(elem_num*9);

  for (int e = 0; e < elem_num; e++){

	const double sqrt_v = sqrt(tetmesh->getElementVolume(e));
	Sqrt_V[e] = sqrt_v;
	for (int i = 0; i < 9; ++i){
	  V.insert(e*9+i,e*9+i) = sqrt_v;
	}
  }
}

/*
 * A = | L C^t F^t|
 *     | C  O  O  |
 *     | F  O  O  |
 * 
 */
void RS2Euler::assembleA(){

  const int n = L.rows() + C.rows() + F.rows();
  const int nz = L.nonZeros() +C.nonZeros()*2 + F.nonZeros()*2;

  typedef Eigen::Triplet<double> T;
  vector<T> A_triplet;
  A_triplet.reserve(nz);

  // assign L to A_triplet
  for (int k=0; k<L.outerSize(); ++k){
	for (SparseMatrix<double>::InnerIterator it(L,k); it; ++it){
	  A_triplet.push_back(  T(it.row(),it.col(),it.value())  );
	}
  }

  // assign C and C^t to A_triplet
  const int Ln = L.rows();
  for (int k=0; k<C.outerSize(); ++k){
	for (SparseMatrix<double>::InnerIterator it(C,k); it; ++it){
	  A_triplet.push_back( T(it.row() + Ln,it.col(), it.value()) );
	  A_triplet.push_back( T(it.col(),it.row() + Ln, it.value()) );
	}
  }

  // assign F and F^t to A_triplet
  const int LnCn = Ln+C.rows();
  for (int k=0; k<F.outerSize(); ++k){
	for (SparseMatrix<double>::InnerIterator it(F,k); it; ++it){
	  A_triplet.push_back( T(it.row() + LnCn,it.col(), it.value()) );
	  A_triplet.push_back( T(it.col(),it.row() + LnCn, it.value()) );
	}
  }

  // set A from the tripletlist
  A.resize (n,n);
  A.reserve (nz);
  A.setFromTriplets(A_triplet.begin(), A_triplet.end());
}

void RS2Euler::assemble_b(const VectorXd &y,VectorXd & b)const{

  const int elem_num = tetmesh->numElements();
  assert_eq(y.size(), elem_num*9);
  b.resize(elem_num*9);

  /// compute b from y
  for (int j = 0; j < elem_num; ++j){
    const double *yj = &(y[9*j]);
	double *bj = &(b[9*j]);
	ComputeBj::compute(yj,Sqrt_V[j],bj);
  }
}
