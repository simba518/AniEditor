#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <ComputeBj.h>
#include <RSCoordComp.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(RSCoordCompTest)

BOOST_AUTO_TEST_CASE(constructWithWarpTest){
  
  const int elem_num = 1;
  const int node_num = 4;
  SparseMatrix<double> G(elem_num*9, node_num*3);
  VectorXd p(node_num*3), y(elem_num*9);

  for (int i = 0; i < G.rows(); ++i){
    for (int j =0; j < G.cols(); ++j){
	  G.insert(i,j) = i*10 + j;
	}
  }

  for (int i = 0; i < p.size(); ++i){
    p[i] = i*15;
  }

  RSCoordComp::constructWithWarp(G, p, y);

  const VectorXd v = G*p;
  Matrix<double, 3, 3, RowMajor> m(&v[0]), M;
  
  M(0,0) = y[0+3];  M(1,1) = y[3+3];   M(2,2) = y[5+3]; 
  M(0,1) = y[1+3]-y[2];  M(0,2) = y[2+3]+y[1];   M(1,2) = y[4+3]-y[0]; 
  M(1,0) = y[1+3]+y[2];  M(2,0) = y[2+3]-y[1];   M(2,1) = y[4+3]+y[0]; 
  
  ASSERT_EQ_SMALL_MAT (m, M);
}

BOOST_AUTO_TEST_CASE(constructWithoutWarpTest){
	  
  const int elem_num = 1;

  // init G as an identity matrix for convinience
  SparseMatrix<double> G(elem_num*9, elem_num*9);
  for (int i = 0; i < G.rows(); ++i){
	G.insert(i,i) = 1.0f;
  }

  // no rotation
  VectorXd u(elem_num*9), y(elem_num*9), correct_y(elem_num*9);
  u << 1,0,0,    0,2,0,   0,0,3;
  RSCoordComp::constructWithoutWarp(G, u, y);
  correct_y << 0,0,0,1,0,0,2,0,3;
  ASSERT_EQ_SMALL_VEC(correct_y,y,elem_num*9);

  // with rotation
  u << 1,0.8,0,    0,2,0.5,   0,0,3;
  RSCoordComp::constructWithoutWarp(G, u, y); // m ==> y

  // // y ==> m: m = exp(y_w)(I+y_s) - I
  Matrix<double, 3, 3, RowMajor> m;
  ComputeBj::compute(&y[0], 1.0f, &m(0,0));

  // comparison
  const VectorXd v = G*u;
  const Matrix<double, 3, 3, RowMajor> old_m(&v[0]);
  ASSERT_EQ_SMALL_MAT_TOL (m,old_m,1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
