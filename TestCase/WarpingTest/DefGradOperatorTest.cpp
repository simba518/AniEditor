#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <DefGradOperator.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(DefGradOperatorTest)

BOOST_AUTO_TEST_CASE(assembleGe2GTest){
  
  typedef Eigen::Triplet<double> T;
  const int elem_v[4] = {0,1,3,2};
  const int e = 1;
  const MatrixXd mGe = MatrixXd::Random(9,12);
  double Ge[9][12];
  for (int i = 0; i < 9; ++i){
    for (int j =0; j < 12; ++j){
	  Ge[i][j] = mGe(i,j);
	}
  }
  vector<T> G_T;
  DefGradOperator::assembleGe2G(elem_v,e,Ge,G_T);
  
  SparseMatrix<double> G(2*9,5*3);
  G.setFromTriplets(G_T.begin(),G_T.end());
  
  MatrixXd correct_G(2*9,5*3);
  correct_G.setZero();
  correct_G.block(9, elem_v[0]*3, 9, 3) = mGe.block(0,0, 9, 3);
  correct_G.block(9, elem_v[1]*3, 9, 3) = mGe.block(0,3, 9, 3);
  correct_G.block(9, elem_v[2]*3, 9, 3) = mGe.block(0,6, 9, 3);
  correct_G.block(9, elem_v[3]*3, 9, 3) = mGe.block(0,9, 9, 3);

  const MatrixXd dG = G;
  ASSERT_EQ_SMALL_MAT(dG,correct_G);
}

BOOST_AUTO_TEST_SUITE_END()
