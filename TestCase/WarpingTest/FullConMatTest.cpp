#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <FullConMat.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(FullConMatTest)

BOOST_AUTO_TEST_CASE(testinsert_I3x3){
  
  SparseMatrix<double> C(9,12);
  FullConMat::insert_I3x3(C,0,0);
  FullConMat::insert_I3x3(C,3,2);

  SparseMatrix<double> corect_C(9,12);
  corect_C.insert(0,0) = 1;
  corect_C.insert(1,1) = 1;
  corect_C.insert(2,2) = 1;

  corect_C.insert(6,9) = 1;
  corect_C.insert(7,10) = 1;
  corect_C.insert(8,11) = 1;

  const MatrixXd dC = C;
  const MatrixXd dcorect_C = corect_C;
  ASSERT_EQ_SMALL_MAT(dC,dcorect_C);
}

BOOST_AUTO_TEST_CASE(testCompute){
	  
  vector<int> con_nodes;
  con_nodes.push_back(0);
  con_nodes.push_back(1);
  con_nodes.push_back(3);

  SparseMatrix<double> C;
  int total_node_num = 4;
  FullConMat::compute(con_nodes,C,total_node_num);
  
  SparseMatrix<double> corect_C(9,12);
  corect_C.insert(0,0) = 1;
  corect_C.insert(1,1) = 1;
  corect_C.insert(2,2) = 1;

  corect_C.insert(3,3) = 1;
  corect_C.insert(4,4) = 1;
  corect_C.insert(5,5) = 1;

  corect_C.insert(6,9) = 1;
  corect_C.insert(7,10) = 1;
  corect_C.insert(8,11) = 1;

  const MatrixXd dC = C;
  const MatrixXd dcorect_C = corect_C;
  ASSERT_EQ_SMALL_MAT(dC,dcorect_C);
}

BOOST_AUTO_TEST_SUITE_END()
