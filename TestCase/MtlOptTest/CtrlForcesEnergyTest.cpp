#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <CtrlForcesEnergy.h>
#include <MtlOptSimulator.h>
#include <Timer.h>
#include <SpaceTimeHessianAutoDiff.h>
using namespace Eigen;
using namespace UTILITY;
using namespace MTLOPT;
using namespace IEDS;

namespace MTLOPT_TEST{

  namespace CTRLFORCESENERGY_TEST{
    
	const int n = 1000;
	const int cubTets = 50;
	const int rb = 80;
	const int rw = 80;
	const int rs = 28;
	const double penaltyCon = 10.0f;
	const double penaltyS = 10.0f;

	const int T = 300;
	const int fixHead = 0;
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
	}
  };

};

BOOST_AUTO_TEST_SUITE(CtrlForcesEnergyTest)

BOOST_AUTO_TEST_CASE(CtrlForcesEnergyZTest){

  using namespace MTLOPT_TEST::CTRLFORCESENERGY_TEST;

  pMtlOptDM dm = pMtlOptDM(new MtlOptDM);
  initMtlOptDM(dm);
  
  CtrlForcesEnergyZ energyZ(dm);
  energyZ.initHessStruct();
  energyZ.reset(true);
  const VectorXd subZ = VectorXd::Random(dm->subFrames()*dm->reducedDim());
  const MatrixXd subZM = Map<MatrixXd>(const_cast<double*>(&subZ[0]), dm->reducedDim(), dm->subFrames());

  { // check timings

  	SparseMatrix<double> HLower;
  	Timer timer;
  	timer.start();
  	for (int i = 0; i < 100; ++i){
  	  energyZ.initHessStruct(HLower);
  	}
  	timer.stop("init struct of H: ");

  	timer.start();
  	for (int i = 0; i < 100; ++i){
  	  energyZ.reset(false);
  	  energyZ.hessian(HLower);
  	}
  	timer.stop("set data of H: ");

  	VectorXd g;
  	timer.start();
  	for (int i = 0;  i< 100; ++i){
	  g = HLower.selfadjointView<Lower>()*subZ;
  	}
  	timer.stop("gradient: ");

  	timer.start();
  	for (int i = 0;  i< 100; ++i){
	  dm->updateZ(&subZM(0,0));
  	  const double fun = energyZ.fun();
  	}
  	timer.stop("fun: ");
  }

  { // check values
	VectorXi fixHeads(5);
	VectorXi fixTails(5);
	fixHeads << 0,1,2,3,4;
	fixTails << 0,1,2,3,4;
	for (int j = 0; j < fixTails.size(); ++j){

	  dm->fixTailFrames(fixTails[j]);

	  for (int i = 0; i < fixHeads.size(); ++i){

		dm->fixHeadFrames(fixHeads[i]);

		const VectorXd subZ = VectorXd::Random(dm->subFrames()*dm->reducedDim());
		const MatrixXd subZM = Map<MatrixXd>(const_cast<double*>(&subZ[0]), dm->reducedDim(), dm->subFrames());
		{
		  // hessian
		  energyZ.initHessStruct();
		  energyZ.reset(true);
		  const SparseMatrix<double> H = energyZ.getHessian().selfadjointView<Lower>();
		  SpaceTimeHessianAutoDiff HessianComp;
		  HessianComp.setTotalFrame(dm->T);
		  HessianComp.setTimeStep(1.0f);
		  HessianComp.setMassDamping(dm->alphaM0);
		  HessianComp.setStiffnessDamping(dm->alphaK0);
		  HessianComp.setEigenValues(dm->Lambda);

		  SparseMatrix<double> Had_full;
		  TEST_ASSERT( HessianComp.compute_H(Had_full) );

		  const int r = dm->reducedDim();
		  const int t0 = dm->T_begin;
		  const int t1 = dm->T_end;
		  const int len = (t1-t0+1)*r;
		  const SparseMatrix<double> Had = Had_full.block(t0*r,t0*r,len,len);
		  ASSERT_EQ(H.rows(),Had.rows());
		  ASSERT_EQ(H.cols(),Had.cols());
		  assert_eq(H.nonZeros(),Had.nonZeros()); ///@bug failed when T=300,rs=30.
		  const double hdiff_norm = (H-Had).norm();
		  ASSERT_EQ_TOL(hdiff_norm, 0.0f, 1e-10);

		  // fun
		  dm->updateZ(&(subZM(0,0)));
		  const VectorXd Z = Map<VectorXd>(&(dm->Z(0,0)), dm->Z.size());
		  const double ztHz_ad = HessianComp.funValue(Z);
		  const double fun = energyZ.fun();
		  ASSERT_EQ_TOL( ztHz_ad, fun ,1e-12*fun);

		  // // grad
		  // const VectorXd g = H*subZ;
		  // VectorXd grad;
		  // energyZ.grad(subZ,grad);
		  // ASSERT_EQ( g, grad);
		}
	  }
	}
  }
}

BOOST_AUTO_TEST_CASE(MechanicalEnergyTest){
  
  using namespace MTLOPT_TEST::CTRLFORCESENERGY_TEST;
  pMtlOptDM dm = pMtlOptDM(new MtlOptDM);
  initMtlOptDM(dm);

  const double penalty_v = 2.0f;
  const double penalty_z = 3.0f;
  MechanicalEnergy energy(dm);
  dm->Lambda = VectorXd::Random(dm->reducedDim())+VectorXd::Ones(dm->reducedDim())*10.0f;
  energy.setPenalty(penalty_v, penalty_z);
  energy.reset(2);

  const VectorXd subZ = VectorXd::Random(dm->subFrames()*dm->reducedDim());
  dm->updateZ(&subZ[0]);

  const SparseMatrix<double> &HLower = energy.getHessian();
  const SparseMatrix<double> H = HLower.selfadjointView<Lower>();
  const double f = energy.fun();
  const double ztHz = subZ.dot(H*subZ)*0.5f;
  ASSERT_EQ_TOL( f, ztHz, 1e-12 );
}

BOOST_AUTO_TEST_CASE(MtlOptSimulatorTest){
 
  using namespace MTLOPT_TEST::CTRLFORCESENERGY_TEST;
  pMtlOptDM dm = pMtlOptDM(new MtlOptDM);
  initMtlOptDM(dm);

  MtlOptSimulatorImplicit simulator(dm);
  simulator.setTotalFrames(100);
  simulator.simulate();

}

BOOST_AUTO_TEST_SUITE_END()
