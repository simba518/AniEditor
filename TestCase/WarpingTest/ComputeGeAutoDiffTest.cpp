#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <ComputeGeAutoDiff.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(ComputeGeAutoDiffTest)

BOOST_AUTO_TEST_CASE(computeTest){

  const double V[4][3] = {{-1,-1,-1},   {3,1,1},  {1,2,1},  {1,2,3} };
  SparseMatrix<double> sGe(9,12);

  double Ge[9][12];
  memset(&Ge[0][0],0,9*12*sizeof(double));
  ComputeGeAutoDiff::getInstance()->compute(V,Ge);
  for (int i = 0; i < 9; ++i){
    for (int j = 0; j < 12; ++j){
	  sGe.insert(i,j) = Ge[i][j];
	}
  }

  VectorXd rest_v(12);
  rest_v << -1,-1,-1,   3,1,1,  1,2,1,  1,2,3;
  for (int i = 0; i < rest_v.size(); ++i){
    rest_v[i] += 10.0f;
  }

  VectorXd I(9);
  I << 1,0,0,  0,1,0,  0,0,1;
  const VectorXd m = sGe*rest_v;
  ASSERT_EQ_SMALL_VEC(I, m, 9);
}

BOOST_AUTO_TEST_SUITE_END()
