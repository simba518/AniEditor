#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <PosConEnergy.h>
#include <MtlOptEnergyAD.h>
#include <CtrlForcesEnergy.h>
#include <CtrlForcesEnergyAD.h>
#include <BasisSelRegEnergy.h>
#include <Timer.h>
#include <CASADITools.h>
using namespace Eigen;
using namespace UTILITY;
using namespace MTLOPT;

namespace MTLOPT_TEST{
  
  const int n = 10000;
  const int cubTets = 100;
  const int rb = 80;
  const int rw = 80;
  const int rs = 3;
  const double penaltyCon = 10.0f;
  const double penaltyS = 10.0f;

  const int T = 20;
  const int fixHead = 2;
  const int fixTail = 2;
  const double h = 1.0f;
  const double ak = 0.3f;
  const double am = 0.1f;
  const VectorXd initialLambda = (VectorXd::Ones(rw)*10+VectorXd::Random(rw));
  const VectorXd lambda = (VectorXd::Ones(rs)*10+VectorXd::Random(rs))*10.0f;
  
};

void initialMtlOptDM_for_test(pMtlOptDM dm){

  using namespace MTLOPT_TEST;

  dm->setTotalFrames(T);
  dm->setDamping(ak,am);
  dm->fixHeadFrames(fixHead);
  dm->fixTailFrames(fixTail);

  dm->Z = MatrixXd::Random(rs,T);
  dm->Lambda = lambda;
  dm->initVariables(lambda,rs);
  dm->S = MatrixXd::Random(rw,rs);

  dm->Z = MatrixXd::Random(rs,T);
  dm->Z.leftCols(fixHead).setZero();
  dm->Z.rightCols(fixTail).setZero();

  dm->conFrames.clear();
  dm->conNodes.clear();
  dm->uc.clear();
  VectorXi conF;
  conF.resize(10);
  for (int f = 0; f < conF.size(); ++f){
    conF[f] = fixHead+f;
  }
  for (int i = 0; i < conF.size(); ++i){

    dm->conFrames.push_back(conF[i]);
  	vector<int> nodes;
	for (int j = 0; j < 100; ++j){
	  nodes.push_back(i+j);
	}
  	dm->conNodes.push_back(nodes);
  	dm->uc.push_back(VectorXd::Random(nodes.size()*3));
  }
}

void initialRedRS_for_test(pReducedRS rs){

  using namespace MTLOPT_TEST;

  const MatrixXd B = MatrixXd::Random(n*3,rb);
  const MatrixXd P = MatrixXd::Random(rb,cubTets*9);
  const MatrixXd hatW = MatrixXd::Random(cubTets*9,rw);
  const VectorXd vol = VectorXd::Ones(cubTets*9)*10.0f+VectorXd::Random(cubTets*9);
  const MatrixXd y = MatrixXd::Random(cubTets*9,T)*0.0f;

  rs->setBP(B,P);
  rs->setHatW(hatW);
  rs->setSqrtV(vol,true);
  rs->setRefY(y);
  rs->checkDimensions();
}

void initialRedSpaceTimeEnergyAD_for_test(RedSpaceTimeEnergyAD &ad){
  
  using namespace MTLOPT_TEST;
  ad.setT(T);
  ad.setTimestep(h);
  ad.setDamping(ak,am,lambda);
  ad.setK(lambda);

  vector<int> fixedFrames;
  for (int i = 0; i < fixHead; ++i){
    fixedFrames.push_back(i);
  }
  for (int i = 0; i < fixTail; ++i){
    fixedFrames.push_back(T-i-1);
  }
  MatrixXd Fz(rs,fixedFrames.size());
  Fz.setZero();
  ad.setFixframes(Fz,fixedFrames);
  
  ad.usePartialCon(true);
  ad.setPartialConPenalty(penaltyCon);
  ad.setFullConPenalty(0.0f);
}

BOOST_AUTO_TEST_SUITE(PosConEnergyTest)

BOOST_AUTO_TEST_CASE(testPosConZ){

  pReducedRS rs = pReducedRS(new ReducedRS());
  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  
  initialRedRS_for_test(rs);
  initialMtlOptDM_for_test(dm);
  
  PosConEnergyZ poscon_z(dm,rs);
  poscon_z.setPenalty(MTLOPT_TEST::penaltyCon);
  poscon_z.updateConstraints();
  poscon_z.reset();
  
  double fun;
  const MatrixXd subZ = MatrixXd::Random(dm->reducedDim(),dm->subFrames());

  VectorXd grad(subZ.size());
  grad.setZero();

  // check time
  {
  	Timer timer; 
  	timer.start();
  	for (int i = 0; i < 200; ++i){
  	  poscon_z.updateYc(&subZ(0,0));
  	  poscon_z.funAddGrad(&fun,&grad[0]);
  	}
  	timer.stop();
  }

  // check value
  {
  	grad.setZero();
  	poscon_z.updateYc(&subZ(0,0));
  	poscon_z.funAddGrad(&fun,&grad[0]);

  	pRedRSWarperAD redad = pRedRSWarperAD(new RedRSWarperAD());
  	redad->setBP(rs->getB(),rs->getP());
  	redad->setHatW(rs->getHatW());
  	redad->setSqrtV(rs->getSqrtV());
  	redad->checkDimensions();

  	RedSpaceTimeEnergyAD ad;
  	initialRedSpaceTimeEnergyAD_for_test(ad);
  	ad.setWSelectMatrix(dm->S);
  	ad.setInitialLambda(MTLOPT_TEST::initialLambda);
  	ad.setPartialCon(dm->conFrames, dm->conNodes, dm->uc);
  	ad.setSPenalty(0.0f);
  	ad.setCtrlForcesPenalty(0.0f);
  	ad.setWarper(redad);
  	ad.assembleEnergy();

  	// check function value
  	const SXMatrix &energy = ad.getEnergy();
  	CasADi::SXFunction funAd = CasADi::SXFunction(ad.getVarZ(),energy);
  	funAd.init();
  	funAd.setInput(&(subZ(0,0)));
  	funAd.evaluate();
  	const CasADi::DMatrix out = funAd.output();
  	ASSERT_EQ_TOL ( out.elem(0,0),fun, (1e-5)*fun );

  	// check gradient
  	const CasADi::SXMatrix gradMat = funAd.jac();
  	CasADi::SXFunction gradFun = CasADi::SXFunction(ad.getVarZ(),gradMat);
  	gradFun.init();
  	VectorXd grad_ad;
  	const VectorXd &subZv = Map<VectorXd>(const_cast<double*>(&subZ(0,0)),subZ.size());
  	CASADI::evaluate(gradFun,subZv,grad_ad);

  	ASSERT_EQ(grad_ad.size(), grad.size());
  	ASSERT_EQ_SMALL_VEC_TOL((grad_ad.transpose()), (grad.transpose()),grad.size(),(1e-5)*grad.norm());

  	// check Jacbian
  	MatrixXd J;
  	poscon_z.updateYc(&subZ(0,0));
  	poscon_z.jacbian(J);

  	VectorXd jp(subZ.size());
  	jp.setZero();
  	int p = 0;
  	for (int i = 0; i < dm->conNodes.size(); ++i){

  	  const int f = dm->conFrames[i] - dm->T_begin;
  	  VectorXd u;
  	  redad->warp(dm->S*subZ.col(f), u, dm->conNodes[i]);
  	  ASSERT_EQ(u.size(), dm->uc[i].size());
  	  const MatrixXd &Jit = J.block(p,0,u.size(),J.cols()).transpose();
  	  jp.segment(f*J.cols(),J.cols()) = Jit*(u-dm->uc[i])*sqrt(MTLOPT_TEST::penaltyCon);
  	  p += u.size();
  	}
  	ASSERT_EQ(jp.size(), grad_ad.size());
  	ASSERT_EQ_SMALL_VEC_TOL((grad_ad.transpose()), (jp.transpose()),jp.size(),(1e-5)*jp.norm());
  }

  // check JtJfull structure
  {
	const int rs = 2;
	const int Ts = 5;
	vector<int> conf;
	conf.push_back(2);
	conf.push_back(4);

	SparseMatrix<double> JtJ;
	poscon_z.initJtJfullStruct(JtJ, conf, rs, Ts, 1);
	ASSERT_EQ( JtJ.nonZeros(), 8 );
	for (int i = 0; i < JtJ.nonZeros(); ++i){
	  JtJ.valuePtr()[i] = i+1;
	}

	vector< Eigen::Triplet<double>  > v;
	v.push_back(Eigen::Triplet<double>(2,2,1));
	v.push_back(Eigen::Triplet<double>(3,3,4));
	v.push_back(Eigen::Triplet<double>(3,2,2));
	v.push_back(Eigen::Triplet<double>(2,3,3));
	
	v.push_back(Eigen::Triplet<double>(6,6,5));
	v.push_back(Eigen::Triplet<double>(7,7,8));
	v.push_back(Eigen::Triplet<double>(7,6,6));
	v.push_back(Eigen::Triplet<double>(6,7,7));

	SparseMatrix<double> JtJ_correct;
	JtJ_correct.resize(10,10);
	JtJ_correct.setFromTriplets(v.begin(), v.end());
	
	ASSERT_EQ((JtJ_correct-JtJ).norm(), 0.0f);
  }

  // check JtJfull
  {
	SparseMatrix<double> JtJ;
	MatrixXd J;
	poscon_z.updateYc(&subZ(0,0));
	poscon_z.jacbian(J);
	poscon_z.initJtJfullStruct(JtJ);
	poscon_z.assembleJtJfull(J, JtJ);
	const SparseMatrix<double> JtJt = JtJ.transpose();
	ASSERT_EQ( (JtJt-JtJ).norm(), 0.0f );
  }
}

BOOST_AUTO_TEST_CASE(testPosConS){

  pReducedRS rs = pReducedRS(new ReducedRS());
  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  
  initialRedRS_for_test(rs);
  initialMtlOptDM_for_test(dm);
  
  PosConEnergyS poscon_s(dm,rs);
  poscon_s.setPenalty(MTLOPT_TEST::penaltyCon);
  poscon_s.updateConstraints();
  poscon_s.reset();
  
  double fun;
  const MatrixXd &S = dm->S;

  VectorXd grad(S.size());
  grad.setZero();

  // check time
  {
  	Timer timer; 
  	timer.start();
  	for (int i = 0; i < 200; ++i){
	  poscon_s.updateYc(&S(0,0));
  	  poscon_s.funAddGrad(&fun,&grad[0]);
  	}
  	timer.stop();
  }

  // check value
  {
  	grad.setZero();
  	poscon_s.updateYc(&S(0,0));
  	poscon_s.funAddGrad(&fun,&grad[0]);

  	pRedRSWarperAD redad = pRedRSWarperAD(new RedRSWarperAD());
  	redad->setBP(rs->getB(),rs->getP());
  	redad->setHatW(rs->getHatW());
  	redad->setSqrtV(rs->getSqrtV());
  	redad->checkDimensions();

  	RedSpaceTimeEnergyAD ad;
  	initialRedSpaceTimeEnergyAD_for_test(ad);
	{ // fix all frames
	  vector<int> fixedFrames;
	  for (int i = 0; i < dm->T; ++i){
		fixedFrames.push_back(i);
	  }
	  ad.setFixframes(dm->Z,fixedFrames);
	}
	
	const VSX sx = makeSymbolic(dm->S.size(),"s");
  	ad.setWSelectMatrix( CASADI::convert(sx,dm->S.cols()) );

  	ad.setInitialLambda(MTLOPT_TEST::initialLambda);
  	ad.setPartialCon(dm->conFrames, dm->conNodes, dm->uc);
  	ad.setSPenalty(0);
  	ad.setCtrlForcesPenalty(0.0f);
  	ad.setWarper(redad);
  	ad.assembleEnergy();

  	// check function value
  	const SXMatrix &energy = ad.getEnergy();
  	CasADi::SXFunction funAd = CasADi::SXFunction(sx,energy);
  	funAd.init();
  	funAd.setInput(&(S(0,0)));
  	funAd.evaluate();
  	const CasADi::DMatrix out = funAd.output();
  	ASSERT_EQ_TOL ( out.elem(0,0),fun, (1e-8)*fun );

  	// check gradient
  	const CasADi::SXMatrix gradMat = funAd.jac();
  	CasADi::SXFunction gradFun = CasADi::SXFunction(sx,gradMat);
  	gradFun.init();
  	VectorXd grad_ad;
  	const VectorXd &subZv = Map<VectorXd>(const_cast<double*>(&S(0,0)),S.size());
  	CASADI::evaluate(gradFun,subZv,grad_ad);

  	ASSERT_EQ(grad_ad.size(), grad.size());
  	ASSERT_EQ_SMALL_VEC_TOL((grad_ad.transpose()), (grad.transpose()),grad.size(),(1e-5)*grad.norm());

  	// check Jacbian
  	MatrixXd J;
  	poscon_s.updateYc(&S(0,0));
  	poscon_s.jacbian(J);

  	VectorXd jp(S.size());
  	jp.setZero();
  	int p = 0;
  	for (int i = 0; i < dm->conNodes.size(); ++i){

  	  const int f = dm->conFrames[i];
  	  VectorXd u;
  	  redad->warp(dm->S*dm->Z.col(f), u, dm->conNodes[i]);
  	  ASSERT_EQ(u.size(), dm->uc[i].size());
  	  const MatrixXd &Jit = J.block(p,0,u.size(),J.cols()).transpose();
  	  jp += Jit*(u-dm->uc[i])*sqrt(MTLOPT_TEST::penaltyCon);
  	  p += u.size();
  	}

  	ASSERT_EQ(jp.size(), grad_ad.size());
  	ASSERT_EQ_SMALL_VEC_TOL((grad_ad.transpose()),(jp.transpose()),jp.size(),(1e-8)*jp.norm());
  }
}

BOOST_AUTO_TEST_CASE(testCtrlForcesEnergyLambdaDamping){

  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  initialMtlOptDM_for_test(dm);

  CtrlForcesEnergyLambdaDamping lambdaEngy(dm);
  lambdaEngy.reset();

  // check value
  {
	const int r = dm->reducedDim();
	const VectorXd lambda_d = VectorXd::Random(r*2);

  	RedSpaceTimeEnergyAD ad;
  	initialRedSpaceTimeEnergyAD_for_test(ad);
  	{ // fix all frames
  	  vector<int> fixedFrames;
  	  for (int i = 0; i < dm->T; ++i){
  		fixedFrames.push_back(i);
  	  }
  	  ad.setFixframes(dm->Z,fixedFrames);
  	}
	
  	const VSX lambdaX = makeSymbolic(r,"lambda");
  	const VSX dampingX = makeSymbolic(r,"d");
	VSX lambdaX_dampingX;
	lambdaX_dampingX.insert(lambdaX_dampingX.end(), lambdaX.begin(), lambdaX.end());
	lambdaX_dampingX.insert(lambdaX_dampingX.end(), dampingX.begin(), dampingX.end());

  	ad.setK(lambdaX);
	ad.setDamping(dampingX);
  	ad.setSPenalty(0.0f);
  	ad.setPartialConPenalty(0.0f);
  	ad.setCtrlForcesPenalty(1.0f);
  	ad.assembleEnergy();

  	// check function value
	const double fun = lambdaEngy.fun(&lambda_d[0]);

  	const SXMatrix &energy = ad.getEnergy();
  	CasADi::SXFunction funAd = CasADi::SXFunction(lambdaX_dampingX,energy);
  	funAd.init();
  	funAd.setInput(&(lambda_d[0]));
  	funAd.evaluate();
  	const CasADi::DMatrix out = funAd.output();
  	ASSERT_EQ_TOL ( out.elem(0,0),fun, (1e-5)*fun );

  	// check gradient
	VectorXd grad(r*2);
	lambdaEngy.grad(&lambda_d[0], &grad[0]);

  	const CasADi::SXMatrix gradMat = funAd.jac();
  	CasADi::SXFunction gradFun = CasADi::SXFunction(lambdaX_dampingX,gradMat);
  	gradFun.init();
  	VectorXd grad_ad;
  	CASADI::evaluate(gradFun,lambda_d,grad_ad);

  	ASSERT_EQ(grad_ad.size(), grad.size());
  	ASSERT_EQ_SMALL_VEC_TOL((grad_ad.transpose()), (grad.transpose()),grad.size(),(1e-5)*grad.norm());

	// check hessian
	const CasADi::SXMatrix hessMat = funAd.hess();
  	CasADi::SXFunction hessFun = CasADi::SXFunction(lambdaX_dampingX,hessMat);
  	hessFun.init();
  	MatrixXd hess_ad;
  	CASADI::evaluate(hessFun,lambda_d,hess_ad);
	
	SparseMatrix<double> hess;
	lambdaEngy.hessianLower(hess);
	const SparseMatrix<double> hs = hess.selfadjointView<Lower>();
	const MatrixXd H = hs;
	const double diff = (hess_ad - H).norm();
	ASSERT_LT(diff,1e-14);

	vector<int> iRow;
	vector<int> iCol;
	vector<double> value;
	lambdaEngy.hessianLower(iRow, iCol, value);
	MatrixXd He;
	He.resize(H.rows(),H.cols());
	He.setZero();
	for (int i = 0; i< value.size(); ++i){
	  He(iRow[i], iCol[i]) = value[i];
	}
	He -= hess;
	ASSERT_EQ(He.norm(),0.0f);

	VectorXd H_lambda, H_damping;
	lambdaEngy.hessianLambda(H_lambda);
	lambdaEngy.hessianDamping(H_damping);
	
	// cout << hess << endl;
	// cout << H_lambda.transpose() << endl;
	// cout << H_damping.transpose() << endl;

  }
}

BOOST_AUTO_TEST_CASE(testBasisSelRegEnergy){

  pReducedRS rs = pReducedRS(new ReducedRS());
  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  
  initialRedRS_for_test(rs);
  initialMtlOptDM_for_test(dm);

  BasisSelRegEnergy basis_sel;
  basis_sel.setPenalty(MTLOPT_TEST::penaltyS);

  double fun;
  const MatrixXd &S = dm->S;

  VectorXd grad(S.size());
  grad.setZero();
  
  // check value
  {
  	pRedRSWarperAD redad = pRedRSWarperAD(new RedRSWarperAD());
  	redad->setBP(rs->getB(),rs->getP());
  	redad->setHatW(rs->getHatW());
  	redad->setSqrtV(rs->getSqrtV());
  	redad->checkDimensions();

  	RedSpaceTimeEnergyAD ad;
  	initialRedSpaceTimeEnergyAD_for_test(ad);
  	{ // fix all frames
  	  vector<int> fixedFrames;
  	  for (int i = 0; i < dm->T; ++i){
  		fixedFrames.push_back(i);
  	  }
  	  ad.setFixframes(dm->Z,fixedFrames);
  	}
	
  	const VSX sx = makeSymbolic(dm->S.size(),"s");
  	ad.setWSelectMatrix( CASADI::convert(sx,dm->S.cols()) );

  	ad.setInitialLambda(VectorXd::Ones(MTLOPT_TEST::initialLambda.size()));
  	ad.setPartialCon(dm->conFrames, dm->conNodes, dm->uc);
  	ad.setPartialConPenalty(0.0f);
	ad.setSPenalty(MTLOPT_TEST::penaltyS);
  	ad.setCtrlForcesPenalty(0.0f);
  	ad.setWarper(redad);
  	ad.assembleEnergy();

  	// check function value
  	const SXMatrix &energy = ad.getEnergy();
  	CasADi::SXFunction funAd = CasADi::SXFunction(sx,energy);
  	funAd.init();
  	funAd.setInput(&(S(0,0)));
  	funAd.evaluate();
  	const CasADi::DMatrix out = funAd.output();

	fun = basis_sel.fun(&S(0,0), S.rows(), S.cols());
  	ASSERT_EQ_TOL ( out.elem(0,0),fun, (1e-10)*fun );

  	// check gradient
  	const CasADi::SXMatrix gradMat = funAd.jac();
  	CasADi::SXFunction gradFun = CasADi::SXFunction(sx,gradMat);
  	gradFun.init();
  	VectorXd grad_ad;
  	const VectorXd &subZv = Map<VectorXd>(const_cast<double*>(&S(0,0)),S.size());
  	CASADI::evaluate(gradFun,subZv,grad_ad);

  	ASSERT_EQ(grad_ad.size(), grad.size());
	basis_sel.grad( &S(0,0), &grad[0], S.rows(), S.cols());
  	ASSERT_EQ_SMALL_VEC_TOL((grad_ad.transpose()), (grad.transpose()),grad.size(),(1e-10)*grad.norm());

	// check hessian
	const CasADi::SXMatrix hessMat = funAd.hess();
  	CasADi::SXFunction hessFun = CasADi::SXFunction(sx,hessMat);
  	hessFun.init();
  	MatrixXd hess_ad;
  	CASADI::evaluate(hessFun,subZv,hess_ad);
	
	VectorXd hess;
	basis_sel.hessian(hess, S.rows(), S.cols());
	const MatrixXd H = hess.asDiagonal();
	const double diff = (hess_ad - H).norm();
	ASSERT_LT(diff,1e-12);
  }
}

BOOST_AUTO_TEST_CASE(testCtrlForcesEnergyAtA_ad){

  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  initialMtlOptDM_for_test(dm);

  CtrlForcesEnergyAtA_AD AtAEngy(dm);
  AtAEngy.reset();

  // check value
  {
	const int r = dm->reducedDim();
	const MatrixXd A = MatrixXd::Random(r,r);
  	const double fun = AtAEngy.fun(&A(0,0), A.size());

  	VectorXd grad(r*r);
  	AtAEngy.grad(&A(0,0), A.size(), &grad[0]);

	MatrixXd H(A.size(), A.size());
	AtAEngy.hessian(&A(0,0), A.size(), &H(0,0));

  	RedSpaceTimeEnergyAD ad;
  	initialRedSpaceTimeEnergyAD_for_test(ad);
  	{ // fix all frames
  	  vector<int> fixedFrames;
  	  for (int i = 0; i < dm->T; ++i){
  		fixedFrames.push_back(i);
  	  }
  	  ad.setFixframes(dm->Z,fixedFrames);
  	}
	
  	const VSX Ax = makeSymbolic(A.size(),"a");
	const SXMatrix As = convert(Ax,r);
	const SXMatrix Ks = CasADi::trans(As).mul(As);
  	ad.setK(Ks);
	SXMatrix D = Ks*dm->alphaK0;
	for (int i = 0; i < D.size1(); ++i){
	  D.elem(i,i) += dm->alphaM0;
	}
  	ad.setDamping(D);

  	ad.setSPenalty(0.0f);
  	ad.setPartialConPenalty(0.0f);
  	ad.setCtrlForcesPenalty(1.0f);
  	ad.assembleEnergy();

  	// check function value
  	const SXMatrix &energy = ad.getEnergy();
  	CasADi::SXFunction funAd = CasADi::SXFunction(Ax,energy);
  	funAd.init();
  	funAd.setInput(&A(0,0));
  	funAd.evaluate();
  	const CasADi::DMatrix out = funAd.output();
  	ASSERT_EQ_TOL ( out.elem(0,0),fun, 1e-12 );

  	// check gradient
  	const CasADi::SXMatrix gradMat = funAd.jac();
  	CasADi::SXFunction gradFun = CasADi::SXFunction(Ax,gradMat);
  	gradFun.init();
  	VectorXd grad_ad;
	const VectorXd Av = Map<VectorXd>(const_cast<double*>(&A(0,0)),A.size());
  	CASADI::evaluate(gradFun,Av,grad_ad);
  	ASSERT_EQ(grad_ad.size(), grad.size());
	ASSERT_LE((grad_ad-grad).norm(), 1e-12);

  	// check hessian
  	const CasADi::SXMatrix hessMat = funAd.hess();
  	CasADi::SXFunction hessFun = CasADi::SXFunction(Ax,hessMat);
  	hessFun.init();
  	MatrixXd hess_ad;
  	CASADI::evaluate(hessFun,Av,hess_ad);
	ASSERT_LE((hess_ad-H).norm(), 1e-12);

  }
}

BOOST_AUTO_TEST_CASE(testCtrlForcesEnergyLambda){

  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  initialMtlOptDM_for_test(dm);

  CtrlForcesEnergyLambda lambdaEngy(dm);
  lambdaEngy.reset();

  // check value
  {
	const VectorXd &lambda = dm->Lambda;
	const double fun = lambdaEngy.fun(lambda);
	VectorXd grad_0;
  	lambdaEngy.grad(grad_0);

  	RedSpaceTimeEnergyAD ad;
  	initialRedSpaceTimeEnergyAD_for_test(ad);
  	{ // fix all frames
  	  vector<int> fixedFrames;
  	  for (int i = 0; i < dm->T; ++i){
  		fixedFrames.push_back(i);
  	  }
  	  ad.setFixframes(dm->Z,fixedFrames);
  	}
	
  	const VSX lambdaX = makeSymbolic(lambda.size(),"lambda");
  	ad.setK(lambdaX);
	ad.setDamping(dm->alphaK0,dm->alphaM0,lambdaX);
  	ad.setSPenalty(0.0f);
  	ad.setPartialConPenalty(0.0f);
  	ad.setCtrlForcesPenalty(1.0f);
  	ad.assembleEnergy();

  	// check function value
  	const SXMatrix &energy = ad.getEnergy();
  	CasADi::SXFunction funAd = CasADi::SXFunction(lambdaX,energy);
  	funAd.init();
  	funAd.setInput(&(lambda[0]));
  	funAd.evaluate();
  	const CasADi::DMatrix out = funAd.output();
  	ASSERT_EQ_TOL ( out.elem(0,0),fun, (1e-5)*fun );

  	// check gradient
  	const CasADi::SXMatrix gradMat = funAd.jac();
  	CasADi::SXFunction gradFun = CasADi::SXFunction(lambdaX,gradMat);
  	gradFun.init();
  	VectorXd grad_ad;
	VectorXd lambda0(lambda.size());
	lambda0.setZero();
  	CASADI::evaluate(gradFun,lambda0,grad_ad);

  	ASSERT_EQ(grad_ad.size(), grad_0.size());
  	ASSERT_EQ_SMALL_VEC_TOL((grad_ad.transpose()), (grad_0.transpose()),grad_0.size(),(1e-5)*grad_0.norm());

	// check hessian
	const CasADi::SXMatrix hessMat = funAd.hess();
  	CasADi::SXFunction hessFun = CasADi::SXFunction(lambdaX,hessMat);
  	hessFun.init();
  	MatrixXd hess_ad;
  	CASADI::evaluate(hessFun,lambda0,hess_ad);
	
	VectorXd hess;
	lambdaEngy.hessian(hess);
	const MatrixXd H = hess.asDiagonal();
	const double diff = (hess_ad - H).norm();
	ASSERT_LT(diff,1e-9);
  }
}

BOOST_AUTO_TEST_CASE(testCtrlForcesEnergyDamping){
  
  pMtlOptDM dm = pMtlOptDM(new MtlOptDM());
  initialMtlOptDM_for_test(dm);

  CtrlForcesEnergyDamping dampingEngy(dm);
  dampingEngy.reset();

  // check value
  {

	ASSERT_GT(dm->Z.norm(), 0.0f);
	ASSERT_GT(dm->Lambda.norm(), 0.0f);

  	RedSpaceTimeEnergyAD ad;
  	initialRedSpaceTimeEnergyAD_for_test(ad);
  	{ // fix all frames
  	  vector<int> fixedFrames;
  	  for (int i = 0; i < dm->T; ++i){
  		fixedFrames.push_back(i);
  	  }
  	  ad.setFixframes(dm->Z,fixedFrames);
  	}

	const VectorXd &lambda = dm->Lambda;
	const int r = lambda.size();
	const VSX Am = makeSymbolic(r,"am");
	const VSX Ak = makeSymbolic(r,"ak");
	VSX AmAk;
	AmAk.insert(AmAk.end(), Am.begin(), Am.end());
	AmAk.insert(AmAk.end(), Ak.begin(), Ak.end());

  	ad.setK(lambda);
	ad.setDamping(Ak,Am,lambda);
  	ad.setSPenalty(0.0f);
  	ad.setPartialConPenalty(0.0f);
  	ad.setCtrlForcesPenalty(1.0f);
  	ad.assembleEnergy();

  	// check function value
	VectorXd am_ak = VectorXd::Random(r*2)*10.0f;
	const VectorXd am = am_ak.head(r);
	const VectorXd ak = am_ak.tail(r);
	const double fun = dampingEngy.fun(am, ak);

  	const SXMatrix &energy = ad.getEnergy();
  	CasADi::SXFunction funAd = CasADi::SXFunction(AmAk,energy);
  	funAd.init();
  	funAd.setInput(&(am_ak[0]));
  	funAd.evaluate();
  	const CasADi::DMatrix out = funAd.output();
  	ASSERT_EQ_TOL ( out.elem(0,0),fun, (1e-8)*fun );

  	// check gradient
	VectorXd grad_0;
  	dampingEngy.grad(grad_0);

  	const CasADi::SXMatrix gradMat = funAd.jac();
  	CasADi::SXFunction gradFun = CasADi::SXFunction(AmAk,gradMat);
  	gradFun.init();
  	VectorXd grad_ad;
	am_ak.setZero();
  	CASADI::evaluate(gradFun,am_ak,grad_ad);

  	ASSERT_EQ(grad_ad.size(), grad_0.size());
  	ASSERT_EQ_SMALL_VEC_TOL((grad_ad.transpose()), (grad_0.transpose()),grad_0.size(),(1e-5)*grad_0.norm()); /// @bug not passed, reason is unknown.

	// check hessian
	MatrixXd hess;
	dampingEngy.hessian(hess);

	const CasADi::SXMatrix hessMat = funAd.hess();
  	CasADi::SXFunction hessFun = CasADi::SXFunction(AmAk,hessMat);
  	hessFun.init();
  	MatrixXd hess_ad;
  	CASADI::evaluate(hessFun,am_ak,hess_ad);

	// cout << hess_ad << endl;
	// cout << hess << endl;

	const double diff = (hess_ad - hess).norm();
	ASSERT_LT(diff,1e-9);
  }
}

BOOST_AUTO_TEST_SUITE_END()
