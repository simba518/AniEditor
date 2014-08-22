#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <KeyframeConEnergy.h>
#include <CASADITools.h>
using namespace Eigen;
using namespace CASADI;
using namespace MTLOPT;

namespace MTLOPT_TEST{

  namespace KEYFRAMECONENERGY_TEST{
    
	const int n = 1000;
	const int cubTets = 50;
	const int rb = 80;
	const int rw = 80;
	const int rs = 2;
	const double penaltyCon = 10.0f;
	const double penaltyS = 10.0f;

	const int T = 7;
	const int fixHead = 2;
	const int fixTail = 0;
	const double h = 1.0f;
	const double ak = 0.01;
	const double am = 0.03;
	const VectorXd initialLambda = VectorXd::Ones(rw)*10+VectorXd::Random(rw);
	const VectorXd lambda = (VectorXd::Ones(rs)*10+VectorXd::Random(rs));

	void initMtlOptDM(pMtlOptDM dm){

	  dm->setTotalFrames(T);
	  dm->setDamping(ak,am);
	  dm->fixHeadFrames(fixHead);
	  dm->fixTailFrames(fixTail);

	  dm->Z = MatrixXd::Random(rs,T);
	  dm->Z.leftCols(fixHead).setZero();
	  dm->Z.rightCols(fixTail).setZero();

	  dm->initVariables(lambda,rs);
	  dm->S = MatrixXd::Random(rw,rs);

	  dm->conFrames.clear();
	  dm->conNodes.clear();
	  dm->uc.clear();
	  VectorXi conF;
	  conF.resize(2);
	  for (int f = 0; f < conF.size(); ++f){
		conF[f] = fixHead+f;
	  }
	  for (int i = 0; i < conF.size(); ++i){

		dm->conFrames.push_back(conF[i]);
		vector<int> nodes;
		for (int j = 0; j < 2; ++j){
		  nodes.push_back(i+j);
		}
		dm->conNodes.push_back(nodes);
		dm->uc.push_back(VectorXd::Random(nodes.size()*3)*0.0f);
	  }
	  
	  vector<int> keyf;
	  keyf.push_back(3);
	  keyf.push_back(5);
	  const MatrixXd Zk = MatrixXd::Random(dm->S.rows(),keyf.size());
	  dm->setKeyfames(keyf,Zk);
	  dm->Z = MatrixXd::Random(rs,dm->T);
	}
  };

};

BOOST_AUTO_TEST_SUITE(KeyframeConEnergyTest)

BOOST_AUTO_TEST_CASE(testEnergyS){

  using namespace MTLOPT_TEST::KEYFRAMECONENERGY_TEST;
  const double penalty = 3.0f;
  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  initMtlOptDM(dm);
  KeyframeConEnergyS ene(dm);
  ene.setPenalty(penalty);
  ene.reset();
  const double fun = ene.fun(&(dm->S(0,0)));
  ASSERT_GT(fun,0.0f);

  const VSX s = makeSymbolic(dm->S.size(),"s");
  const CasADi::SXMatrix sm = convert(s, dm->S.cols());
  const CasADi::SXMatrix zk = convert(dm->Zk);
  CasADi::SXMatrix zm(dm->Z.rows(),dm->Zk.cols());

  for (int i = 0; i < dm->Zk.cols(); ++i){
	for (int j = 0; j < dm->Z.rows(); ++j){
	  zm.elem(j,i) = dm->Z(j,dm->keyframes[i]);
	}
  }
  ASSERT_EQ(zm.size1(),sm.size2());
  ASSERT_EQ(sm.size1(),zk.size1());
  const CasADi::SXMatrix sz_zk = sm.mul(zm)-zk;

  CasADi::SXMatrix E = 0.0f;
  for (int i = 0; i < sz_zk.size1(); ++i){
	for (int j = 0; j < sz_zk.size2(); ++j)
	  E += sz_zk.elem(i,j)*sz_zk.elem(i,j)*penalty;
  }
  E *= 0.5f;

  // check value
  const VectorXd sx = Map<VectorXd>(const_cast<double*>(&(dm->S(0,0))),dm->S.size());
  CasADi::SXFunction funAd = CasADi::SXFunction(s,E);
  funAd.init();
  funAd.setInput(&(sx[0]));
  funAd.evaluate();
  const CasADi::DMatrix out = funAd.output();
  ASSERT_EQ_TOL ( out.elem(0,0),fun, fun*1e-12);
  
  // check grad
  VectorXd g(sx.size());
  g.setZero();
  ene.gradAdd(&(sx[0]),&(g[0]));

  const CasADi::SXMatrix gradMat = funAd.jac();
  CasADi::SXFunction gradFun = CasADi::SXFunction(s,gradMat);
  gradFun.init();
  VectorXd grad_ad;
  CASADI::evaluate(gradFun,sx,grad_ad);
  ASSERT_LT((grad_ad-g).norm(),1e-12);
}

BOOST_AUTO_TEST_CASE(testEnergyZ){

  using namespace MTLOPT_TEST::KEYFRAMECONENERGY_TEST;
  const double penalty = 3.0f;
  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  initMtlOptDM(dm);
  KeyframeConEnergyZ ene(dm);
  ene.setPenalty(penalty);
  ene.updateConstraints();
  const double fun = ene.fun(&(dm->Z(0,dm->T_begin)));
  ASSERT_GT(fun,0.0f);

  const VSX z = makeSymbolic(dm->subFrames()*dm->reducedDim(),"z");
  const CasADi::SXMatrix zk = convert(dm->Zk);
  const CasADi::SXMatrix sm = convert(dm->S);
  CasADi::SXMatrix zm(sm.size2(),zk.size2());
  
  for (int i = 0; i < zm.size1(); ++i){
  	for (int j = 0; j < zm.size2(); ++j){
	  ASSERT_GT(dm->keyframes[j]-dm->T_begin,0);
	  ASSERT_GT(z.size(),(dm->keyframes[j]-dm->T_begin)*dm->reducedDim()+i);
  	  zm.elem(i,j) = z[(dm->keyframes[j]-dm->T_begin)*dm->reducedDim()+i];
	}
  }

  ASSERT_EQ(zm.size1(),sm.size2());
  ASSERT_EQ(sm.size1(),zk.size1());
  const CasADi::SXMatrix sz_zk = sm.mul(zm)-zk;

  CasADi::SXMatrix E = 0.0f;
  for (int i = 0; i < sz_zk.size1(); ++i){
  	for (int j = 0; j < sz_zk.size2(); ++j)
  	  E += sz_zk.elem(i,j)*sz_zk.elem(i,j)*penalty;
  }
  E *= 0.5f;

  // check value
  const VectorXd zx = Map<VectorXd>(const_cast<double*>(&(dm->Z(0,dm->T_begin))),z.size());
  CasADi::SXFunction funAd = CasADi::SXFunction(z,E);
  funAd.init();
  funAd.setInput(&(zx[0]));
  funAd.evaluate();
  const CasADi::DMatrix out = funAd.output();
  ASSERT_EQ_TOL ( out.elem(0,0),fun, fun*1e-12);
  
  // check grad
  VectorXd g(zx.size());
  g.setZero();
  ene.gradAdd(&(zx[0]),&(g[0]));

  const CasADi::SXMatrix gradMat = funAd.jac();
  CasADi::SXFunction gradFun = CasADi::SXFunction(z,gradMat);
  gradFun.init();
  VectorXd grad_ad;
  CASADI::evaluate(gradFun,zx,grad_ad);
  ASSERT_LT((grad_ad-g).norm(),1e-12);
  
  // check hessian
  const CasADi::SXMatrix hesMat = funAd.hess();
  CasADi::SXFunction hesFun = CasADi::SXFunction(z,hesMat);
  hesFun.init();

  MatrixXd hes_ad;
  CASADI::evaluate(hesFun,zx,hes_ad);
  const MatrixXd hh = ene.getHess();
  const MatrixXd d  = hh-hes_ad;
  ASSERT_LT(d.norm(),1e-12);

}

BOOST_AUTO_TEST_SUITE_END()
