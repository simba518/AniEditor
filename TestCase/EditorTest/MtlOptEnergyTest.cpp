#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MtlOptEnergy.h>
#include <MtlOptIpopt.h>
#include <CASADITools.h>
#include <MatrixTools.h>
#include <MASimulatorAD.h>
#include <Timer.h>
#include <JsonFilePaser.h>
using namespace Eigen;
using namespace UTILITY;
using namespace LSW_ANI_EDITOR;

class TestIpoptObjfun:public BaseFunGrad{
  // (x0-2)^2 + (x1-1)^2
public:
  int dim()const{
	return 2;
  }
  void init(double *x,const int n){
	x[0] = 0;
	x[1] = 0;
  }
  double fun(const double *x){
	return (x[0]-2)*(x[0]-2) + (x[1]-1)*(x[1]-1);
  }
  void grad(const double *x,double *g){
	g[0] = 2*(x[0]-2);
	g[1] = 2*(x[1]-1);
  }
  void setRlst(const double *x, const double objValue){
	_objval = objValue;
	_x[0] = x[0];
	_x[1] = x[1];
  }

  const Vector2d &getOptX()const{
	return _x;
  }
  double getObjVal()const{
	return _objval;
  }

private:
  double _objval;
  Vector2d _x;
};
typedef boost::shared_ptr<TestIpoptObjfun> pTestIpoptObjfun;

typedef NoWarp TemptTestWaper;
typedef boost::shared_ptr<TemptTestWaper> pTemptTestWaper;
typedef CtrlForceEnergyT<pTemptTestWaper> CtrlForceEnergyTempt;
typedef boost::shared_ptr<CtrlForceEnergyTempt> pCtrlForceEnergyTempt;

class MtlOptEnergyTestDataModel{
  
public:
  MtlOptEnergyTestDataModel(){

	{// input data
	  h = 0.3;
	  r = 2;
	  T = 4;
	  alpha_k = 0.1;
	  alpha_m = 0.2;	  

	  lambda = VectorXd::Random(r)+VectorXd::Ones(r)*10.0f;
	  v0 = VectorXd::Random(r)*1.2f;
	  z0 = VectorXd::Random(r);
	  w = VectorXd::Random(r*T-r);
	}

	{// constraints
	  penaltyCon = 100.0f;
	  penaltyKey = 1000.0f;
	  conF.push_back(1);
	  conF.push_back(3);
	  conN.push_back(vector<int>(1,1));
	  conN.push_back(vector<int>(1,2));
	  conN[conN.size()-1].push_back(3);
	  for (int i = 0; i < conN.size(); ++i)
		uc.push_back(VectorXd::Random(conN[i].size()*3));

	  keyframes.push_back(1);
	  keyframes.push_back(3);
	  for (int i = 0; i < keyframes.size(); ++i){
		keyZ.push_back(VectorXd::Random(lambda.size()));
	  }
	}
	
	{ // set all as zero
	  // v0.setZero();
	  // z0.setZero();
	  // w.setZero();
	  // for (int i = 0; i < uc.size(); ++i)
	  // 	uc[i].setZero();
	}

	{// initialize
	  warper = pTemptTestWaper(new TemptTestWaper(MatrixXd::Random(30*3,r)));
	  ctrlF = pCtrlForceEnergyTempt(new CtrlForceEnergyTempt());
	  ctrlF->setRedWarper(warper);

	  ctrlF->setTimestep(h);
	  ctrlF->setTotalFrames(T);
	  const VectorXd diagD = alpha_k*lambda+VectorXd::Ones(r)*alpha_m;
	  ctrlF->setMtl(lambda,diagD);
	  ctrlF->setV0(v0);
	  ctrlF->setZ0(z0);
	  ctrlF->precompute();

	  ctrlF->setPenaltyCon(penaltyCon,penaltyKey);
	  ctrlF->setPartialCon(conF,conN,uc);
	  ctrlF->setKeyframes(keyZ,keyframes);
	}
  }
  void initSimulator(MASimulatorAD &simulator){
	
	simulator.setTimeStep(h);
	simulator.setEigenValues<VectorXd>(lambda);
	simulator.setStiffnessDamping(SX(alpha_k));
	simulator.setMassDamping(SX(alpha_m));
	simulator.setIntialStatus(v0,z0);
  }
  void initMtlEnergy(MtlOptEnergy &mtl){
	mtl.setTotalFrames(T);
	mtl.setTimestep(h);
	mtl.setMtl(lambda,alpha_k,alpha_m);
  }
  
public:
  int r;
  int T;
  double h;
  double alpha_k;
  double alpha_m;
  double penaltyCon;
  double penaltyKey;

  VectorXd lambda,v0,z0,w;

  vector<int> keyframes;
  vector<VectorXd> keyZ;
  vector<int> conF;
  vector<vector<int> > conN;
  vector<VectorXd> uc;

  pTemptTestWaper warper;
  pCtrlForceEnergyTempt ctrlF;
};

template<typename WARPER_POINTER>
class CtrlForceEnergyTestADT{
  
public:
  void setRedWarper(WARPER_POINTER warper){
	_warper = warper;
  }
  MASimulatorAD &getSimulator(){
	return _simulator;
  }
  void setTotalFrames(const int T){
	_T = T;
  }
  void setPenaltyCon(const double pc,const double pk){
	assert_ge(pc,0.0f);
	assert_ge(pk,0.0f);
	_penaltyCon = pc;
	_penaltyKey = pk;
  }
  void setKeyframes(const vector<VectorXd> &kZ, const vector<int> &keyf){
	_keyZ = kZ;
	_keyframes = keyf;
  }
  void setPartialCon(vector<int>&conF,vector<vector<int> >&conN,vector<VectorXd>&uc){
	_conFrames = conF;
	_conNodes = conN;
	_uc = uc;
  }
  
  void precompute(){

	_w = CASADI::makeSymbolic(dim(),"w");
	SXMatrix E = 0.0f;
	for (int i = 0; i < _w.size(); ++i){
	  E += _w[i]*_w[i];
	}
	const vector<SXMatrix> w = CASADI::Sx2Mat(_w,_T-1);
	vector<SXMatrix> V,Z;
	_simulator.forward(w,V,Z);

	{ // full constraints
	  for (size_t i = 0; i < _keyframes.size(); ++i){
		const SXMatrix &zi = Z[_keyframes[i]];
		const SXMatrix zk = CASADI::convert((VectorXd)_keyZ[i]);
		assert_eq(zk.size1(),zi.size1());
		assert_eq(zk.size2(),zi.size2());
		assert_eq(zk.size2(),1);
		for (int k = 0; k < zi.size(); ++k){
		  const SX n = (zi.elem(k,0)-zk.elem(k,0));
		  E += _penaltyKey*n*n;
		}		
	  }
	}

	{ // partial constraints
	  for (int i = 0; i < _conFrames.size(); ++i){
		const int f = _conFrames[i];
		const SXMatrix ui = _warper->warp(Z[f],_conNodes[i]);
		const SXMatrix suc = CASADI::convert((VectorXd)_uc[i]);
		assert_eq(ui.size1(),suc.size1());
		assert_eq(ui.size2(),suc.size2());
		for (int k = 0; k < ui.size(); ++k){
		  const SX n = (ui.elem(k,0)-suc.elem(k,0));
		  E += _penaltyCon*n*n;
		}
	  }
	}

	E = 0.5f*E;
	_energyFun = CasADi::SXFunction(_w,E);
	_energyFun.init();
	const SXMatrix g = _energyFun.jac();
	_gradFun = CasADi::SXFunction(_w,g);
	_gradFun.init();
  }

  double fun(const double *x){
	_energyFun.setInput(x);
	_energyFun.evaluate();
	double f = 0.0f;
	_energyFun.getOutput(&f);
	return f;
  }
  void grad(const double *x,double *g){
	_gradFun.setInput(x);
	_gradFun.evaluate();
	_gradFun.getOutput(g);
  }
  int dim()const{
	return (_T-1)*reducedDim();
  }
  int reducedDim()const{
	return _simulator.reducedDim();
  }
  
private:
  int _T;
  MASimulatorAD _simulator;
  WARPER_POINTER _warper;

  double _penaltyCon;
  double _penaltyKey;
  vector<int> _keyframes;
  vector<VectorXd> _keyZ;
  vector<int> _conFrames;
  vector<vector<int> > _conNodes;
  vector<VectorXd> _uc;

  CasADi::SXFunction _energyFun;
  CasADi::SXFunction _gradFun;
  vector<SX> _w;
};
typedef CtrlForceEnergyTestADT<pTemptTestWaper> CtrlForceEnergyTestAD;

class MtlOptEnergyAD:public MtlOptEnergy{
  
public:
  void precompute(){

	const int r = reducedDim();
	const vector<SX> akama = CASADI::makeSymbolic(dim(),"x");
	const SXMatrix ak = CASADI::makeEyeMatrix(vector<SX>(akama.begin(),akama.begin()+r));
	const SXMatrix am = CASADI::makeEyeMatrix(vector<SX>(akama.begin()+r,akama.begin()+r*2));
	const SXMatrix A = CASADI::convert(vector<SX>(akama.begin()+2*r,akama.end()),r);
	const SXMatrix K = CasADi::trans(A).mul(A);
	const SXMatrix D = am+ak.mul(K);
	
	SXMatrix E = 0.0f;
	const int T = _Z.size();
	for (int k = 0; k < T-1; ++k){
	  const SXMatrix vk = CASADI::convert(_V[k]);
	  const SXMatrix vk1 = CASADI::convert(_V[k+1]);
	  const SXMatrix zk1 = CASADI::convert(_Z[k+1]);
	  const SXMatrix n = vk1-vk+_h*D.mul(vk1)+_h*K.mul(zk1);
	  for (int i = 0; i < r; ++i)
		E += n.elem(i,0)*n.elem(i,0);
	}

	E = 0.5f*E;
	_energyFun = CasADi::SXFunction(akama,E);
	_energyFun.init();
	const SXMatrix g = _energyFun.jac();
	_gradFun = CasADi::SXFunction(akama,g);
	_gradFun.init();
  }
  double fun(const double *x){
	_energyFun.setInput(x);
	_energyFun.evaluate();
	double f = 0.0f;
	_energyFun.getOutput(&f);
	return f;
  }
  void grad(const double *x,double *g){
	_gradFun.setInput(x);
	_gradFun.evaluate();
	_gradFun.getOutput(g);
  }
  void setVZ(const MatrixXd &V,const MatrixXd &Z){
	MtlOptEnergy::setVZ(V,Z);
	EIGEN3EXT::convert(V,_V);
	EIGEN3EXT::convert(Z,_Z);
  }
  
protected:
  
private:
  CasADi::SXFunction _energyFun;
  CasADi::SXFunction _gradFun;
  vector<VectorXd> _V;
  vector<VectorXd> _Z;
};

BOOST_AUTO_TEST_SUITE(MtlOptEnergyTest)

BOOST_AUTO_TEST_CASE(testTemptTestWarper){
 
  TemptTestWaper warper(MatrixXd::Random(30,4));
  const VectorXd z = VectorXd::Random(4);
  vector<int> nodes(2);
  nodes[0] = 1;
  nodes[1] = 3;
  VectorXd ui;
  warper.warp(z,0,nodes,ui);
  
  const SXMatrix sz = CASADI::convert(z);
  const SXMatrix sui = warper.warp(sz,nodes);
  const VectorXd ssui = CASADI::convert2Vec<double>(sui);

  ASSERT_EQ(ssui.size(),ui.size());
  ASSERT_EQ_SMALL_VEC_TOL(ssui,ui,ui.size(),1e-12);
}

BOOST_AUTO_TEST_CASE(testSimpleFunOpt){

  pTestIpoptObjfun energy = pTestIpoptObjfun(new TestIpoptObjfun());
  NoConIpoptSolver solver(energy);
  solver.setTol(1e-8);
  solver.setMaxIt(100);
  TEST_ASSERT( solver.initialize() );
  TEST_ASSERT( solver.solve() );
  
  const double optObjVal = 0.0f;
  Vector2d optX;
  optX << 2.0f,1.0f;
  ASSERT_EQ_TOL( energy->getObjVal(),optObjVal,1e-8);
  ASSERT_EQ_SMALL_VEC_TOL( energy->getOptX(),optX,2,1e-8);
}

BOOST_AUTO_TEST_CASE(testCtrlForward){

  MtlOptEnergyTestDataModel data;
 
  MatrixXd V,Z;
  data.ctrlF->forward(data.w,V,Z);
  
  MASimulatorAD simulator;
  data.initSimulator(simulator);

  for (int i = 0; i < data.r; ++i){

  	SX sG[2][2], sS[2];
  	simulator.getImpIntegMat(sG,sS,i);
  	Eigen::Matrix<double,2,2> G;
  	Eigen::Matrix<double,2,1> S;
  	G(0,0) = sG[0][0].getValue();
  	G(0,1) = sG[0][1].getValue();
  	G(1,0) = sG[1][0].getValue();
  	G(1,1) = sG[1][1].getValue();
  	S(0,0) = sS[0].getValue();
  	S(1,0) = sS[1].getValue();

	const MatrixXd G2 = (data.ctrlF->getG().block<2,2>(i*2,0));
	const MatrixXd S2 = (data.ctrlF->getS().block<2,1>(i*2,0));
  	ASSERT_EQ_SMALL_MAT_TOL(G,G2,1e-12);
  	ASSERT_EQ_SMALL_MAT_TOL(S,S2,1e-12);
  }

  vector<VectorXd> w2,v2,z2;
  EIGEN3EXT::convert(Map<MatrixXd>(&data.w[0],data.r,data.T-1),w2);
  simulator.forward(w2,v2,z2);

  MatrixXd V2,Z2;
  EIGEN3EXT::convert(v2,V2);
  EIGEN3EXT::convert(z2,Z2);
  
  ASSERT_EQ(V.size(),V2.size());
  ASSERT_EQ(Z.size(),Z2.size());
  ASSERT_EQ_SMALL_MAT_TOL(V,V2,1e-10);
  ASSERT_EQ_SMALL_MAT_TOL(Z,Z2,1e-10);
    
}

BOOST_AUTO_TEST_CASE(testCtrlObjValue){

  MtlOptEnergyTestDataModel data;
  CtrlForceEnergyTestAD ctrlFAD;
  data.initSimulator(ctrlFAD.getSimulator());
  ctrlFAD.setPenaltyCon(data.penaltyCon,data.penaltyKey);
  ctrlFAD.setPartialCon(data.conF,data.conN,data.uc);
  ctrlFAD.setKeyframes(data.keyZ,data.keyframes);
  ctrlFAD.setTotalFrames(data.T);
  ctrlFAD.setRedWarper(data.warper);
  ctrlFAD.precompute();

  ASSERT_GT (data.w.size(),0);
  const double obj_1 = data.ctrlF->fun(&data.w[0]);
  const double obj_2 = ctrlFAD.fun(&data.w[0]);
  ASSERT_EQ_TOL(obj_1,obj_2,1e-12);

}

BOOST_AUTO_TEST_CASE(testCtrlGrad){

  MtlOptEnergyTestDataModel data;
  CtrlForceEnergyTestAD ctrlFAD;
  data.initSimulator(ctrlFAD.getSimulator());
  ctrlFAD.setPenaltyCon(data.penaltyCon,data.penaltyKey);
  ctrlFAD.setPartialCon(data.conF,data.conN,data.uc);
  ctrlFAD.setKeyframes(data.keyZ,data.keyframes);
  ctrlFAD.setTotalFrames(data.T);
  ctrlFAD.setRedWarper(data.warper);
  ctrlFAD.precompute();

  ASSERT_GT (data.w.size(),0);

  VectorXd g1(ctrlFAD.dim());
  VectorXd g2(ctrlFAD.dim());

  data.ctrlF->grad(&data.w[0],&g1[0]);
  ctrlFAD.grad(&data.w[0],&g2[0]);

  ASSERT_EQ_SMALL_VEC_TOL(g1,g2,g1.size(),1e-10);

}

BOOST_AUTO_TEST_CASE(testCtrlOpt){

  MtlOptEnergyTestDataModel data;
  NoConIpoptSolver solver(data.ctrlF);
  solver.setTol(1e-8);
  solver.setMaxIt(100);
  // solver.setPrintLevel(5);
  TEST_ASSERT( solver.initialize() );
  TEST_ASSERT( solver.solve() );
}

BOOST_AUTO_TEST_CASE(testCtrlOptTiming){

  const int T = 400;
  JsonFilePaser inf;
  { // open json file
	const string infname = "/home/simba/Workspace/AnimationEditing/Data/beam/aniedit.ini";
	TEST_ASSERT( inf.open(infname) );
  }

  pRedRSWarperExt nodeWarper;
  { // init warper
	pRSWarperExt warper = pRSWarperExt(new RSWarperExt());

	string vol_filename;
	TEST_ASSERT( inf.readFilePath("vol_filename",vol_filename) );
	pTetMesh tet_mesh = pTetMesh(new TetMesh());
	TEST_ASSERT( tet_mesh->load(vol_filename) );
	warper->setTetMesh(tet_mesh);

	vector<int> fixed_nodes;
	TEST_ASSERT(inf.readVecFile("fixed_nodes_RS2Euler",fixed_nodes,UTILITY::TEXT));
	warper->setFixedNodes(fixed_nodes);

	MatrixXd W,B;
	TEST_ASSERT( inf.readMatFile("eigen_vectors",W) );
	TEST_ASSERT( inf.readMatFile("nonlinear_basis",B) );
	const MatrixXd Uref(W.rows(),T);
	vector<VectorXd> u_ref;
	EIGEN3EXT::convert(Uref,u_ref);
	warper->setRefU(u_ref);
	TEST_ASSERT( warper->precompute() );

	vector<int> cubP;
	vector<double> cubW;
	TEST_ASSERT( inf.readVecFile("cubaturePoints",cubP,UTILITY::TEXT) );
	TEST_ASSERT( inf.readVecFile("cubatureWeights",cubW) );
	nodeWarper = pRedRSWarperExt(new RedRSWarperExt(warper,B,W,cubP,cubW));
  }
  
  pCtrlForceEnergy ctrlF = pCtrlForceEnergy(new CtrlForceEnergy);
  { // init ctrlF
	double h, alpha_k, alpha_m, penaltyCon,penaltyKey;
	VectorXd lambda;
	TEST_ASSERT( inf.read("h",h));
	TEST_ASSERT( inf.read("alpha_k",alpha_k));
	TEST_ASSERT( inf.read("alpha_m",alpha_m));
	TEST_ASSERT( inf.read("penaltyCon",penaltyCon));
	TEST_ASSERT( inf.read("penaltyKey",penaltyKey));
	TEST_ASSERT( inf.readVecFile("eigen_values",lambda));

	ctrlF->setTimestep(h);
	ctrlF->setTotalFrames(T);
	ctrlF->setPenaltyCon(penaltyCon,penaltyKey);
	const VectorXd D = alpha_m*VectorXd::Ones(lambda.size())+lambda*alpha_k;
	ctrlF->setMtl(lambda,D);
	ctrlF->precompute();
	ctrlF->setRedWarper(nodeWarper);
	ctrlF->setZ0(VectorXd::Random(lambda.size()));
	ctrlF->setV0(VectorXd::Random(lambda.size()));
  }

  { // set constraints
	vector<int> conF;
	vector<vector<int> > conN;
	vector<VectorXd> uc;
	ASSERT_GE(T,200);
	conF.push_back(20);
	conF.push_back(50);
	conF.push_back(100);
	conF.push_back(120);
	conF.push_back(160);
	conF.push_back(190);

	for (int i = 0; i < conF.size(); ++i){
	  vector<int> nodes;
	  for (int j = 0; j < 10; ++j)
		nodes.push_back(j+i);
	  conN.push_back(nodes);
	}

	for (int i = 0; i < conN.size(); ++i)
	  uc.push_back(VectorXd::Random(conN[i].size()*3));

	ctrlF->setPartialCon(conF,conN,uc);
  }

  { // test timing
	const VectorXd w = VectorXd::Random(ctrlF->dim());
	Timer timer;
	timer.start();
	for (int i = 0; i < 1000; ++i)
	  ctrlF->fun(&w[0]);
	timer.stop("compute function value: ");

	VectorXd g(w.size());
	timer.start();
	for (int i = 0; i < 1000; ++i)
	  ctrlF->grad(&w[0],&g[0]);
	timer.stop("compute gradient: ");
  }
}

BOOST_AUTO_TEST_CASE(testMtlEnergyObjValue){
  
  MtlOptEnergyTestDataModel data;
  const MatrixXd V = MatrixXd::Random(data.ctrlF->reducedDim(),data.T);
  const MatrixXd Z = MatrixXd::Random(data.ctrlF->reducedDim(),data.T);

  MtlOptEnergy mtl;
  data.initMtlEnergy(mtl);
  mtl.setVZ(V,Z);
  const VectorXd akamA = VectorXd::Random(mtl.dim());
  const double objV1 = mtl.fun(&akamA[0]);

  MtlOptEnergyAD mtlAD;
  data.initMtlEnergy(mtlAD);
  mtlAD.setVZ(V,Z);
  mtlAD.precompute();
  const double objV2 = mtlAD.fun(&akamA[0]);
  
  ASSERT_EQ_TOL (objV1,objV2,1e-14);
}

BOOST_AUTO_TEST_CASE(testMtlEnergyObjValue2){
  
  MtlOptEnergyTestDataModel data;
  data.ctrlF->clearPartialCon();
  const VectorXd w = VectorXd::Random(data.ctrlF->dim());
  const double objV1 = data.ctrlF->fun(&w[0]);

  MatrixXd V,Z;
  data.ctrlF->forward(w,V,Z);

  MtlOptEnergy mtl;
  data.initMtlEnergy(mtl);
  mtl.setVZ(V,Z);
  const VectorXd akamA = mtl.getX();
  const double objV2 = mtl.fun(&akamA[0])/(data.h*data.h);
  
  ASSERT_EQ_TOL (objV1,objV2,1e-14);
}

BOOST_AUTO_TEST_CASE(testMtlEnergyGrad){
  
  MtlOptEnergyTestDataModel data;
  const MatrixXd V = MatrixXd::Random(data.ctrlF->reducedDim(),data.T);
  const MatrixXd Z = MatrixXd::Random(data.ctrlF->reducedDim(),data.T);

  MtlOptEnergy mtl;
  data.initMtlEnergy(mtl);
  mtl.setVZ(V,Z);
  const VectorXd akamA = VectorXd::Random(mtl.dim());
  VectorXd g1(akamA.size());
  mtl.grad(&akamA[0],&g1[0]);

  MtlOptEnergyAD mtlAD;
  data.initMtlEnergy(mtlAD);
  mtlAD.setVZ(V,Z);
  mtlAD.precompute();
  VectorXd g2(akamA.size());
  mtlAD.grad(&akamA[0],&g2[0]);
  
  ASSERT_EQ_SMALL_VEC_TOL (g1,g2,g1.size(),1e-14);
}

BOOST_AUTO_TEST_SUITE_END()
