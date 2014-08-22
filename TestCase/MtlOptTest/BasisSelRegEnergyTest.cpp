#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <BasisSelRegEnergy.h>
#include <CASADITools.h>
using namespace Eigen;
using namespace CASADI;
using namespace MTLOPT;

BOOST_AUTO_TEST_SUITE(BasisSelRegEnergyTest)

BOOST_AUTO_TEST_CASE(testRegEnergy){

  const double penalty = 3.0f;
  MatrixXd S(3,2);
  S << 1,2,	3,4, 5,6;
  
  BasisSelRegEnergy energy;
  energy.setPenalty(penalty);
  
  const double f = energy.fun(&S(0,0),S.rows(), S.cols());
  ASSERT_EQ_TOL(f,91.0f*0.5*penalty,1e-10);

  MatrixXd g(3,2),Gc(3,2);
  energy.grad(&S(0,0),&g(0,0),S.rows(), S.cols());
  Gc = S*penalty;
  ASSERT_EQ(g,Gc);

  VectorXd H,cH(6);
  energy.hessian(H,S.rows(),S.cols());
  cH << 1,1,1,1,1,1;
  cH *= penalty;
  ASSERT_EQ(H,cH);
}

BOOST_AUTO_TEST_CASE(testBasisNormEnergy){
 
  const double penalty = 2.0f;
  BasisNormalizeEnergy NormLizeS;
  NormLizeS.setPenalty(penalty);

  const MatrixXd S = MatrixXd::Random(4,3);
  const double fun = NormLizeS.fun(&S(0,0), S.rows(), S.cols());
  ASSERT_GT(fun,0.0f);

  VectorXd g(S.size());
  NormLizeS.grad(&S(0,0), &g[0], S.rows(), S.cols());
  ASSERT_GT(g.norm(), 0.0f);

  const VSX s = makeSymbolic(S.size(),"s");
  const CasADi::SXMatrix sm = convert(s, S.cols());
  const CasADi::SXMatrix I = makeEyeMatrix(VectorXd::Ones(S.cols()));
  const CasADi::SXMatrix sts_i = trans(sm).mul(sm)-I;
  
  CasADi::SXMatrix E = 0.0f;
  for (int i = 0; i < sts_i.size1(); ++i){
	for (int j = 0; j < sts_i.size2(); ++j)
	  E += sts_i.elem(i,j)*sts_i.elem(i,j)*penalty;
  }
  E *= 0.5f;

  // check value
  const VectorXd sx = Map<VectorXd>(const_cast<double*>(&S(0,0)),S.size());
  CasADi::SXFunction funAd = CasADi::SXFunction(s,E);
  funAd.init();
  funAd.setInput(&(sx[0]));
  funAd.evaluate();
  const CasADi::DMatrix out = funAd.output();
  ASSERT_EQ_TOL ( out.elem(0,0),fun, fun*1e-12);
  
  // check grad
  const CasADi::SXMatrix gradMat = funAd.jac();
  CasADi::SXFunction gradFun = CasADi::SXFunction(s,gradMat);
  gradFun.init();

  VectorXd grad_ad;
  CASADI::evaluate(gradFun,sx,grad_ad);
  ASSERT_EQ((grad_ad.transpose()), (g.transpose()));
  
}

BOOST_AUTO_TEST_CASE(testBasisNormColEnergy){
 
  const double penalty = 1.0f;
  BasisNormalizeColEnergy NormLizeS;
  NormLizeS.setPenalty(penalty);

  const MatrixXd S = MatrixXd::Random(4,3);
  const double fun = NormLizeS.fun(&S(0,0), S.rows(), S.cols());
  ASSERT_GT(fun,0.0f);

  VectorXd g(S.size());
  NormLizeS.grad(&S(0,0), &g[0], S.rows(), S.cols());
  ASSERT_GT(g.norm(), 0.0f);

  const VSX s = makeSymbolic(S.size(),"s");
  const CasADi::SXMatrix sm = convert(s, S.cols());
  const CasADi::SXMatrix I = makeEyeMatrix(VectorXd::Ones(S.cols()));
  const CasADi::SXMatrix sts_i = trans(sm).mul(sm)-I;
  
  CasADi::SXMatrix E = 0.0f;
  for (int i = 0; i < sts_i.size1(); ++i){
	E += sts_i.elem(i,i)*sts_i.elem(i,i)*penalty;
  }
  E *= 0.5f;

  // check value
  const VectorXd sx = Map<VectorXd>(const_cast<double*>(&S(0,0)),S.size());
  CasADi::SXFunction funAd = CasADi::SXFunction(s,E);
  funAd.init();
  funAd.setInput(&(sx[0]));
  funAd.evaluate();
  const CasADi::DMatrix out = funAd.output();
  ASSERT_EQ_TOL ( out.elem(0,0),fun, fun*1e-12);
  
  // check grad
  const CasADi::SXMatrix gradMat = funAd.jac();
  CasADi::SXFunction gradFun = CasADi::SXFunction(s,gradMat);
  gradFun.init();

  VectorXd grad_ad;
  CASADI::evaluate(gradFun,sx,grad_ad);
  ASSERT_EQ((grad_ad.transpose()), (g.transpose()));
  
}

BOOST_AUTO_TEST_SUITE_END()
