#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <MatrixIO.h>
#include <SpaceTimeHessianAutoDiff.h>
#include <SpaceTimeHessianInverse.h>
using namespace Eigen;
using namespace UTILITY;
using namespace IEDS;

struct SpaceTimeHessianInverseTestInit{

  SpaceTimeHessianInverseTestInit(){

	const int T = 20;
	const double h = 0.03f;
	const double alpha_k = 0.001f;
	const double alpha_m = 0.02f;
	r = 3;
	VectorXd eigen_values(r);
	eigen_values << 1,2,3;
	pBaseSpaceTimeHessian Hcmp = pBaseSpaceTimeHessian(new SpaceTimeHessianAutoDiff());
	Hcmp->setTotalFrame(T);
	Hcmp->setTimeStep(h);
	Hcmp->setStiffnessDamping(alpha_k);
	Hcmp->setMassDamping(alpha_m);
	Hcmp->setEigenValues(eigen_values);
	TEST_ASSERT(Hcmp->compute_H(H));
	ASSERT_EQ(H.rows(),H.cols());
	ASSERT_EQ(H.rows(),r*T);
	H = block(H,2*r,2*r,H.rows()-4*r,H.rows()-4*r);
  }
  SparseMatrix<double> H;
  int r;
};

BOOST_AUTO_TEST_SUITE(SpaceTimeHessianInverseTest)

BOOST_FIXTURE_TEST_CASE(testModeBlock, SpaceTimeHessianInverseTestInit){
  
  vector<MatrixXd> Ms;
  SpaceTimeHessianInverse::getModeBlocks(H,r,Ms);
  SparseMatrix<double> H2;
  SpaceTimeHessianInverse::assembleFromModeBlocks(Ms, H2);
  const MatrixXd Hm = H;
  const MatrixXd H2m = H2;
  ASSERT_EQ_SMALL_MAT(Hm,H2m);
}

BOOST_FIXTURE_TEST_CASE(testInverse, SpaceTimeHessianInverseTestInit){
	  
  SparseMatrix<double> Hinv;
  SpaceTimeHessianInverse::inverse(H,r,Hinv);

  const SparseMatrix<double> sparseHxHinv = (H*Hinv);
  const SparseMatrix<double> I = eye(H.rows(),(double)1.0f);
  ASSERT_EQ_TOL(I.norm(),sparseHxHinv.norm(),1e-8);
  ASSERT_EQ_TOL(I.norm(),sparseHxHinv.norm(),1e-8);
  ASSERT_EQ_TOL((I-sparseHxHinv).norm(),0,1e-8);
}

BOOST_AUTO_TEST_SUITE_END()
