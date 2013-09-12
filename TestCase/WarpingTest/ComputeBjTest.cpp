#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <ComputeBj.h>
using namespace Eigen;
using namespace LSW_WARPING;

BOOST_AUTO_TEST_SUITE(ComputeBjTest)

BOOST_AUTO_TEST_CASE(computeTest){

  VectorXd yj(9), bj(9);
  const double Sqrt_Vj = 0.5f;
  yj << 1,2,3,4,5,6,7,8,9;
  bj.setZero();

  ComputeBj::compute(&yj[0], Sqrt_Vj, &bj[0]);
  Matrix<double, 3, 3, RowMajor> m(&bj[0]), correct_m;

  correct_m << -0.18562,1.47395,1.21579,
	1.56010,1.53761,2.87480,
	4.35514,5.65027,5.84487;
  ASSERT_EQ_SMALL_MAT_TOL(m,correct_m, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
