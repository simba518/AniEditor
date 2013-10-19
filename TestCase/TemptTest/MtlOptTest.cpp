#include "MtlOptModel.h"

BOOST_AUTO_TEST_SUITE(MtlOptTest)

BOOST_AUTO_TEST_CASE(checkRedSpaceTimeEnergyAD){
  
  RedSpaceTimeEnergyAD energy;
  VectorXd lambda(2);
  lambda << 1,2;
  energy.setT(6);
  energy.setTimestep(0.1);
  energy.setDamping(0.2,0.3,lambda);
  energy.setK(lambda);
  VectorXi Kid(2);
  Kid << 0,5;
  MatrixXd Kz(2,2);
  Kz << 1,1,2,2;
  energy.setKeyframes(Kz,Kid);
  energy.assembleEnergy();

  const VSX &z = energy.getVarZ();
  ASSERT_EQ (z.size(),8);
  ASSERT_EQ (energy.reducedDim(),2);
}

BOOST_AUTO_TEST_CASE(Opt_Z){

  MtlOptModel model;
  model.produceSimRlst();
  // model.extrangeKeyframes();

  RedSpaceTimeEnergyAD ad;
  model.initMtlOpt(ad);
  ad.assembleEnergy();
  model.initSolver(ad.getEnergy(),ad.getVarZ());

  const int lenZ = ad.getVarZ().size();
  vector<double> x0(lenZ,0);
  const VectorXd &cZ = Map<VectorXd>(&model.CorrectZ(0,0),model.CorrectZ.size());
  ASSERT_EQ(lenZ,cZ.size());
  for (int i = 0; i < lenZ; ++i)
    x0[i] = cZ[i];
  model.solver.setInput(x0,CasADi::NLP_X_INIT);

  model.solve();
  model.computeEnergy(Map<VectorXd>(&model.CorrectZ(0,0),model.CorrectZ.size()));
  model.saveRlst();
}

BOOST_AUTO_TEST_CASE(Opt_LambdaEq0){

  MtlOptModel model;
  model.produceSimRlst();
  // model.extrangeKeyframes();

  RedSpaceTimeEnergyAD ad;
  model.initMtlOpt(ad);
  ad.setK(VectorXd::Zero(model.redDim()));
  ad.assembleEnergy();
  model.initSolver(ad.getEnergy(),ad.getVarZ());
  model.computeEnergy(Map<VectorXd>(&model.Z(0,0),model.Z.size()));
  model.solve();
  model.saveRlst();
}

BOOST_AUTO_TEST_CASE(Opt_Z_Damping){

  MtlOptModel model;
  model.produceSimRlst();
  model.extrangeKeyframes();

  RedSpaceTimeEnergyAD ad;
  model.initMtlOpt(ad);
  const VSX damping = makeSymbolic(model.redDim(),"d");
  ad.setDamping(damping);

  ad.assembleEnergy();

  VSX varX = ad.getVarZ();
  varX.insert(varX.end(),damping.begin(),damping.end());
  model.initSolver(ad.getEnergy(),varX);

  model.solve();
  model.saveRlst();
}

BOOST_AUTO_TEST_CASE(Opt_Z_Damping_Lambda){

  MtlOptModel model;
  model.produceSimRlst();
  model.extrangeKeyframes();

  RedSpaceTimeEnergyAD ad;
  model.initMtlOpt(ad);

  const VSX damping = makeSymbolic(model.redDim(),"d");
  ad.setDamping(damping);
  const VSX lambda = makeSymbolic(model.redDim(),"La");
  ad.setK(lambda);

  ad.assembleEnergy();
  VSX varX = ad.getVarZ();
  varX.insert(varX.end(), damping.begin(), damping.end());
  varX.insert(varX.end(),lambda.begin(),lambda.end());
  model.initSolver(ad.getEnergy(),varX);

  vector<double> x0(varX.size(),0);
  const int lenZ = ad.getVarZ().size();
  const VectorXd &cZ = Map<VectorXd>(&model.CorrectZ(0,0),model.CorrectZ.size());
  ASSERT_EQ(lenZ,cZ.size());
  for (int i = 0; i < lenZ; ++i)
    x0[i] = cZ[i]*0.8f;
  model.solver.setInput(x0,CasADi::NLP_X_INIT);

  model.solve();
  model.saveRlst();
}

BOOST_AUTO_TEST_CASE(Opt_Z_Damping_K){

  MtlOptModel model;
  model.produceSimRlst();
  model.extrangeKeyframes();

  RedSpaceTimeEnergyAD ad;
  model.initMtlOpt(ad);

  const VSX damping = makeSymbolic(model.redDim(),"d");
  ad.setDamping(damping);

  const int r = model.redDim();
  const VSX ax = makeSymbolic(r*r,"a");
  const SXMatrix A = convert(ax,r);
  const SXMatrix K = CasADi::trans(A).mul(A);
  ad.setK(K);

  ad.assembleEnergy();
  VSX varX = ad.getVarZ();
  varX.insert(varX.end(), damping.begin(), damping.end());
  varX.insert(varX.end(),ax.begin(),ax.end());
  model.initSolver(ad.getEnergy(),varX);

  vector<double> x0(varX.size(),0);
  const int lenZ = ad.getVarZ().size();
  const VectorXd &cZ = Map<VectorXd>(&model.CorrectZ(0,0),model.CorrectZ.size());
  ASSERT_EQ(lenZ,cZ.size());
  for (int i = 0; i < lenZ; ++i)
    x0[i] = cZ[i]*0.8f;
  MatrixXd initA(r,r);
  initA.setZero();
  for (int i = 0; i < r; ++i)
	initA(i,i) = sqrt(model.lambda[i]);
  const VectorXd &vA = Map<VectorXd>(&initA(0,0),initA.size());
  for (int i = lenZ+damping.size(); i < x0.size(); ++i)
	x0[i] = vA[i-lenZ-damping.size()]*0.8f;

  model.solver.setInput(x0,CasADi::NLP_X_INIT);

  vector<double> lowerB(x0.size(),-std::numeric_limits<double>::infinity());
  for (int i = lenZ; i < x0.size()-ax.size(); ++i)
	lowerB[i] = 0.0f;

  model.solver.setInput(lowerB,CasADi::NLP_LBX);

  model.solve();
  model.saveRlst();

  VectorXd rlstAv = model.getOutput().tail(r*r);
  const MatrixXd rlstA = Map<MatrixXd>(&rlstAv[0],r,r);
  const MatrixXd rlstK = rlstA.transpose()*rlstA;
  cout<< "K = " << endl << rlstK << endl;
  EigenSolver<MatrixXd> eigenK(rlstK);
  cout<< "eigen(K): " << eigenK.eigenvalues().transpose() << endl;
}

BOOST_AUTO_TEST_SUITE_END()
