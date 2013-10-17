#include <iostream>
#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <AuxTools.h>
#include <MatrixIO.h>
#include <MASimulatorAD.h>
#include <TetMeshEmbeding.h>
#include <RS2Euler.h>
#include <RSCoordComp.h>
#include <DefGradOperator.h>
#include <MapMA2RS.h>
#include <MeshVtkIO.h>
#include <limits>
#include "MtlOpt.h"
using namespace std;
using namespace EIGEN3EXT;
using namespace UTILITY;
using namespace LSW_WARPING;

typedef struct _MtlOptModel{

  _MtlOptModel(){

	dataDir = "/home/simba/Workspace/AnimationEditor/Data/beam/";
	T = 200;
	h = 0.3f;
	alphaK = 0.001f;
	alphaM = 0.001f;
	Kid.resize(6);
	Kid << 0,12,50,85,140,199;

	testModeId.resize(2);
	testModeId<<0,1;
	z0.resize(testModeId.size());
	z0[0] = 10000;
	z0[1] = 7000;

	loadLambda(dataDir+"eigen_values80.b");
	const VectorXd tlam = lambda;
	lambda.resize(testModeId.size());
	for (int i = 0; i < testModeId.size(); ++i)
	  lambda[i] = tlam[testModeId[i]];
	initVolObj();
  }
  bool loadLambda(const string fname){
	const bool succ = load(fname,lambda);
	assert(succ);
	return succ;
  }
  void initVolObj(){

	const string evect = dataDir+"/eigen_vectors80.b";
	MatrixXd tW;
	TEST_ASSERT( load(evect,tW) );
	W.resize(tW.rows(),testModeId.size());
	for (int i = 0; i < testModeId.size(); ++i){
	  W.col(i) = tW.col(testModeId[i]);
	}	

	TEST_ASSERT( volobj.loadObjMesh(dataDir+"beam.obj"));
	TEST_ASSERT( volobj.loadTetMesh(dataDir+"beam.abq"));
	TEST_ASSERT( volobj.loadWeights(dataDir+"interp-weights.txt"));
	vector<int> fixed_nodes;
	TEST_ASSERT( loadVec(dataDir+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
	rs2euler.setTetMesh(volobj.getTetMesh());
	rs2euler.setFixedNodes(fixed_nodes);
	rs2euler.precompute();
	TEST_ASSERT( DefGradOperator::compute(volobj.getTetMesh(),G) );
	MapMA2RS::computeMapMatPGW(G,W,PGW);	
  }
  void produceSimRlst(){

	const int r = redDim();
	LSW_ANI_EDITOR::MASimulatorAD sim;
	sim.setTimeStep(h);
	sim.setEigenValues(lambda);
	sim.setStiffnessDamping(alphaK);
	sim.setMassDamping(alphaM);
	const VectorXd zero = VectorXd::Zero(r);
	sim.setIntialStatus(zero,z0);

	Z.resize(r,T);
	CorrectZ.resize(r,(T-Kid.size()));
	int ci = 0;
	for (int i = 0; i < T; ++i){
	  Z.col(i) = sim.getEigenZ();
	  if( RedSpaceTimeEnergyAD::isKeyframe(Kid,i) < 0){
		CorrectZ.col(ci) = sim.getEigenZ();
		ci ++;
	  }
	  sim.forward(zero);
	}

	Kz.resize(r,Kid.size());
	for (int i = 0; i < Kid.size(); ++i){
	  const int f = Kid[i];
	  assert_in( f,0,Z.cols() );
	  Kz.col(i) = Z.col(f);
	}
  }
  void extrangeKeyframes(){
	ASSERT_EQ(Kid.size(),6);
	Kid << 0,50,12,85,140,199;
  }
  void initMtlOpt(RedSpaceTimeEnergyAD &ad){
	ad.setT(T);
	ad.setTimestep(h);
	ad.setDamping(alphaK,alphaM,lambda);
	ad.setK(lambda);
	ad.setKeyframes(Kz,Kid);
  }
  void initSolver(const SXMatrix &E, const VSX &x){

	allVars = x;
	fun = CasADi::SXFunction(x,E);
	solver = CasADi::IpoptSolver(fun);
	solver.setOption("generate_hessian",true);
	solver.init();
	
	const int r = redDim();
	const int lenZ = (T-Kid.size())*r;
	vector<double> x0(x.size(),0);
	vector<double> lowerB(x0.size(),-std::numeric_limits<double>::infinity());
	for (int i = lenZ; i < x0.size(); ++i){
	  x0[i] = 0.01f;
	  lowerB[i] = 0.0f;
	}
	solver.setInput(x0,CasADi::NLP_X_INIT);
	solver.setInput(lowerB,CasADi::NLP_LBX);
  }
  void solve(){

	solver.solve();
	double cost;
	solver.getOutput(cost,CasADi::NLP_COST);
	cout<<"input data: " << endl;
	cout<<"lambda: " << lambda.transpose() << endl << endl;
	cout<< endl << "optimal cost: " << cost << endl << endl;
	const VectorXd xx = getOutput();
	const int lenZ = redDim()*(T-Kid.size());
	cout<< "opt mtl: " << xx.tail(xx.size()-lenZ).transpose() << endl << endl;
  }
  VectorXd getOutput(){

	VectorXd x;
	const int r = redDim();
	vector<double> vx(allVars.size());
	solver.getOutput(vx,CasADi::NLP_X_OPT);
	x.resize(vx.size());
	for (size_t i = 0; i < vx.size(); ++i)
	  x[i] = vx[i];
	ASSERT_GE(x.size(),r*(T-Kid.size()));
	return x;
  }
  void getZfromSolver(MatrixXd &Z){

	const VectorXd x = getOutput();
	const int r = redDim();
	Z.resize(r,T);
	int xpos = 0;
	for (size_t i = 0; i < T; ++i){
	  const int k = RedSpaceTimeEnergyAD::isKeyframe(Kid,i);
	  if(k >= 0){
		Z.col(i) = Kz.col(k);
	  }else{
		ASSERT_LE (xpos,x.size()-r);
		Z.col(i) = x.segment(xpos,r);
		xpos = xpos+r;
	  }
	}
  }
  void saveRlst(){
	
	saveMesh(Z,"input");
	saveMesh(Kz,"key");
	MatrixXd newZ;
	getZfromSolver(newZ);
	saveMesh(newZ,"new");
  }
  void saveMesh(const MatrixXd &Z,const string fname){

	for (int i = 0; i < Z.cols(); ++i){
	  const string ffname=dataDir+"tempt/meshes/"+fname+TOSTR(i)+".vtk";
	  saveMeshOneZ(Z.col(i),ffname);
	}
  }
  void saveMeshOneZ(const VectorXd &z,const string fname){

	VectorXd y,u;
	const VectorXd p = W*z;
	RSCoordComp::constructWithWarp(G,p,y);
	TEST_ASSERT ( rs2euler.reconstruct(y,u) );
	volobj.interpolate(u);
	TEST_ASSERT (LSW_SIM::MeshVtkIO::write(volobj.getObjMesh(),fname));
  }
  void computeEnergy(const VectorXd &X){
	VectorXd rlst;
	CASADI::evaluate(fun,X,rlst);
	cout<< "E(X) = " << rlst.transpose() << endl;
  }
  int redDim()const{
	return lambda.size();
  }

  string dataDir;
  MatrixXd W;
  SparseMatrix<double> G;
  MatrixXd PGW;
  TetMeshEmbeding volobj;
  RS2Euler rs2euler;

  int T;
  double h;
  double alphaK;
  double alphaM;
  VectorXd lambda;
  VectorXi testModeId;
  VectorXd z0;
  MatrixXd Z;
  MatrixXd CorrectZ;
  MatrixXd Kz;
  VectorXi Kid;

  CasADi::SXFunction fun;
  CasADi::IpoptSolver solver;
  VSX allVars;

}MtlOptModel;

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

  // const VSX ak = makeSymbolic(model.redDim(),"ak");
  // const VSX am = makeSymbolic(model.redDim(),"am");
  // ad.setDamping(ak,am,model.lambda);

  ad.assembleEnergy();
  VSX varX = ad.getVarZ();

  varX.insert(varX.end(),damping.begin(),damping.end());
  // varX.insert(varX.end(),ak.begin(),ak.end());
  // varX.insert(varX.end(),am.begin(),am.end());

  model.initSolver(ad.getEnergy(),varX);
  model.solve();
  model.saveRlst();
}

BOOST_AUTO_TEST_CASE(Opt_Z_Lambda){

  MtlOptModel model;
  model.produceSimRlst();

  RedSpaceTimeEnergyAD ad;
  model.initMtlOpt(ad);

  const VSX lambda = makeSymbolic(model.redDim(),"La");
  ad.setK(lambda);
  ad.setDamping(0.150809,0.244,model.lambda);

  ad.assembleEnergy();
  VSX varX = ad.getVarZ();
  varX.insert(varX.end(),lambda.begin(),lambda.end());
  model.initSolver(ad.getEnergy(),varX);

  // vector<double> x0(varX.size(),0);
  // const int lenZ = ad.getVarZ().size();
  // const VectorXd &cZ = Map<VectorXd>(&model.CorrectZ(0,0),model.CorrectZ.size());
  // ASSERT_EQ(lenZ,cZ.size());
  // for (int i = 0; i < lenZ; ++i){
  //   x0[i] = cZ[i]*0.8f;
  // }
  // x0[lenZ] = 1.0f;
  // model.solver.setInput(x0,CasADi::NLP_X_INIT);

  model.solve();
  model.saveRlst();
}

BOOST_AUTO_TEST_CASE(Opt_Z_Damping_Lambda){

  MtlOptModel model;
  model.produceSimRlst();
  // model.extrangeKeyframes();

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
  x0[lenZ] = 0.0159731;
  x0[lenZ+1] = 0.0497447;
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
  for (int i = 0; i < r; ++i){
	initA(i,i) = sqrt(model.lambda[i]);
  }
  const VectorXd &vA = Map<VectorXd>(&initA(0,0),initA.size());
  for (int i = lenZ+damping.size(); i < x0.size(); ++i){
	x0[i] = vA[i-lenZ-damping.size()]*0.8f;
  }

  model.solver.setInput(x0,CasADi::NLP_X_INIT);

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
