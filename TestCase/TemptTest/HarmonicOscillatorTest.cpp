#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include "HarmonicOscillator.h"
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(HarmonicOscillatorTest)

BOOST_AUTO_TEST_CASE(generateSequenceTest){

  const double lambda = 0.5;
  const double alpha_k = 0.1;
  const double alpha_m = 0.2;
  const double z0 = 0.5;
  const double v0 = 0.8;

  HarmonicOscillator<double> z(lambda,alpha_k,alpha_m,z0,v0);
  ASSERT_EQ(z(0),0.5);
  ASSERT_EQ_TOL(z(1.83),1.059137778977706,1e-10);
}

BOOST_AUTO_TEST_SUITE_END()
