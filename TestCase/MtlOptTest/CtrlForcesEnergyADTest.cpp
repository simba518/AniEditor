#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <CtrlForcesEnergyAD.h>
#include <CtrlForcesEnergy.h>
using namespace Eigen;
using namespace MTLOPT;

BOOST_AUTO_TEST_SUITE(CtrlForcesEnergyADTest)

// BOOST_AUTO_TEST_CASE(testCtrlForcesK_AtA_AD){
  
//   pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
//   dm->T = 10;
//   dm->alphaM0 = 0.1;
//   dm->alphaK0 = 0.15;
//   const int r = 3;
//   dm->Z = MatrixXd::Random(r, dm->T)*10.0f;

//   const VectorXd lambda_sqrt = VectorXd::Random(dm->reducedDim());
//   VectorXd lambda(lambda_sqrt.size());
//   for (int i = 0; i < lambda.size(); ++i)
//     lambda[i] = lambda_sqrt[i]*lambda_sqrt[i];
//   dm->initVariables(lambda,r);
//   dm->Z = MatrixXd::Random(r, dm->T)*10.0f;

//   CtrlForcesEnergyLambda lambda_fun(dm); 
//   CtrlForcesEnergyK_AD K_fun(dm);
//   CtrlForcesEnergyAtA_AD AtA_fun(dm);

//   lambda_fun.reset();
//   AtA_fun.reset();
//   K_fun.reset();

//   const double fun_lambda = lambda_fun.fun(lambda);
//   ASSERT_GT(fun_lambda, 0);

//   const MatrixXd A = lambda_sqrt.asDiagonal();
//   const double fun_ata = AtA_fun.fun(&A(0,0), A.size());
//   ASSERT_GT(fun_ata, 0);

//   VectorXd K((r+1)*r/2);
//   K.setZero();
//   for (int i = 0; i < r; ++i) K[K_fun.symIndex(i,i)] = lambda[i];
//   const double fun_K = K_fun.fun(&K(0,0),K.size());
//   ASSERT_GT(fun_K, 0);

//   ASSERT_EQ_TOL(fun_ata, fun_lambda, 1e-12);
//   ASSERT_EQ(fun_ata, fun_K);

//   {
// 	CASADI::VSX Xv;
// 	CasADi::SXMatrix K;
// 	const int r = 3;
// 	K_fun.initK(Xv, K, r);
// 	// cout << K << endl;
// 	// for (int i = 0; i < Xv.size(); ++i){
// 	//   cout << Xv[i] << "\t";
// 	// }
// 	// cout << endl;
//   }

  
// }

BOOST_AUTO_TEST_SUITE_END()
