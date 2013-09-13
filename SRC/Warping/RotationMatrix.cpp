#include <SparseMatrixTools.h>
#include "RotationMatrix.h"
#include <ComputeBj.h>
using namespace LSW_WARPING;

void RotationMatrix::compute(const VectorXd &u, SparseMatrix<double>& R)const{

  // compute rs coordinates of each element
  VectorXd RS_u;
  _warper->computeRSCoord(u, RS_u, false);

  // compute the rotational vector of each node.
  const int node_num = u.size()/3;
  vector<Vector3d> w(node_num);
  vector<int> neighbor_ele_of_node(node_num);
  for (int i = 0; i < node_num; ++i){
    w[i].setZero();
	neighbor_ele_of_node[i] = 0;
  }
  for (int ei = 0; ei < _tet_mesh->numElements(); ++ei){
	const MeshElement *ele = _tet_mesh->element(ei);
    for (int j = 0; j < _tet_mesh->numElementVertices(); ++j){
	  const int vi = ele->vertex(j);
	  w[vi][0] += RS_u[ei*9+0];
	  w[vi][1] += RS_u[ei*9+1];
	  w[vi][2] += RS_u[ei*9+2];
	  neighbor_ele_of_node[vi] ++;
	}
  }
  for (int i = 0; i < node_num; ++i){
    w[i] = w[i]/(double)(neighbor_ele_of_node[i]);
  }

  // compute the rotational matrix of each node
  mat3r Ri;
  vector<mat3r> R_vec;
  R_vec.reserve(node_num);
  for (int i = 0; i < node_num; ++i){
	R_vec.push_back(ComputeBj::compute_exp(w[i],Ri));
  }
  EIGEN3EXT::eye(R,R_vec);
}
