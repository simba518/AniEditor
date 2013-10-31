#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MtlOptEnergy.h>
#include <MtlOptIpopt.h>
#include <CASADITools.h>
#include <MatrixTools.h>
#include <MASimulatorAD.h>
using namespace Eigen;
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

class TemptTestWaper{
  
public:
  // u(z) = W*z
  TemptTestWaper(const int nx3,const int r){
	assert_gt(r,0);
	assert_ge(nx3,r);
	assert_eq(nx3%3,0);
	W = MatrixXd::Random(nx3,r);
  }
  void warp(const VectorXd &z,int frame_id,const vector<int> &nodes,VectorXd &ui){
	assert_eq(z.size(),W.cols());
	ui = block(W,nodes)*z;
  }
  void jacobian(const VectorXd &z,int frame_id,int node,Eigen::Matrix<double,3,-1> &J){
	J = W.block(node*3,0,3,W.cols());
  }
  SXMatrix warp(const SXMatrix &z,const vector<int> &nodes){
	assert_eq(z.size1(),W.cols());
	assert_eq(z.size2(),1);
	SXMatrix sWi;  CASADI::convert(this->block(W,nodes),sWi);
	return sWi.mul(z);
  }
  
protected:
  MatrixXd block(const MatrixXd &W,const vector<int> &nodes)const{

	MatrixXd Wi(nodes.size()*3, W.cols());
	for (size_t i = 0; i < nodes.size(); ++i){
	  const int j3 = nodes[i]*3;
	  assert_le(j3,W.rows()-3);
	  Wi.block( i*3,0,3, W.cols() ) = W.block( j3,0,3,W.cols() );
	}
	return Wi;
  }
  
private:
  MatrixXd W;
};
typedef boost::shared_ptr<TemptTestWaper> pTemptTestWaper;
typedef CtrlForceEnergyT<pTemptTestWaper> CtrlForceEnergyTempt;

class MtlOptEnergyTestDataModel{
  
public:
  MtlOptEnergyTestDataModel(){

	{// input data
	  h = 3.0;
	  r = 2;
	  T = 4;
	  alpha_k = 0.1;
	  alpha_m = 0.2;	  

	  lambda.resize(2);
	  v0.resize(2);
	  z0.resize(2);
	  lambda << 0.25f, 2.5f;
	  v0 << 1,2;
	  z0 << 3,4;
	  w = VectorXd::Random(r*T-r);
	}

	{// constraints
	  penaltyCon = 10.0f;
	  // conF.push_back(1);
	  conF.push_back(3);
	  conN.push_back(vector<int>(1,1));
	  // conN.push_back(vector<int>(1,2));
	  uc.push_back(VectorXd::Random(3));
	  // uc.push_back(VectorXd::Random(3));
	}

	{// initialize
	  warper = pTemptTestWaper(new TemptTestWaper(30*3,r));
	  ctrlF.setRedWarper(warper);

	  ctrlF.setTimestep(h);
	  ctrlF.setTotalFrames(T);
	  const VectorXd diagD = alpha_k*lambda+VectorXd::Ones(r)*alpha_m;
	  ctrlF.setMtl(lambda,diagD);
	  ctrlF.setV0(v0);
	  ctrlF.setZ0(z0);
	  ctrlF.precompute();

	  ctrlF.setPenaltyCon(penaltyCon);
	  ctrlF.setPartialCon(conF,conN,uc);
	}
  }
  void initSimulator(MASimulatorAD &simulator){
	
	simulator.setTimeStep(h);
	simulator.setEigenValues<VectorXd>(lambda);
	simulator.setStiffnessDamping(SX(alpha_k));
	simulator.setMassDamping(SX(alpha_m));
	simulator.setIntialStatus(v0,z0);
  }
  
public:
  int r;
  int T;
  double h;
  double alpha_k;
  double alpha_m;
  double penaltyCon;

  VectorXd lambda,v0,z0,w;

  vector<int> conF;
  vector<vector<int> > conN;
  vector<VectorXd> uc;

  pTemptTestWaper warper;
  CtrlForceEnergyTempt ctrlF;
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
  void setPenaltyCon(const double pc){
	assert_ge(pc,0.0f);
	_penaltyCon = pc;
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
  vector<int> _conFrames;
  vector<vector<int> > _conNodes;
  vector<VectorXd> _uc;

  CasADi::SXFunction _energyFun;
  CasADi::SXFunction _gradFun;
  vector<SX> _w;
};
typedef CtrlForceEnergyTestADT<pTemptTestWaper> CtrlForceEnergyTestAD;


BOOST_AUTO_TEST_SUITE(MtlOptEnergyTest)

BOOST_AUTO_TEST_CASE(testTemptTestWarper){
 
  TemptTestWaper warper(30,4);
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
  data.ctrlF.forward(data.w,V,Z);
  
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

  	ASSERT_EQ(G,(data.ctrlF.getG().block<2,2>(i*2,0)));
  	ASSERT_EQ(S,(data.ctrlF.getS().block<2,1>(i*2,0)));
  }

  vector<VectorXd> w2,v2,z2;
  EIGEN3EXT::convert(Map<MatrixXd>(&data.w[0],data.r,data.T-1),w2);
  simulator.forward(w2,v2,z2);

  MatrixXd V2,Z2;
  EIGEN3EXT::convert(v2,V2);
  EIGEN3EXT::convert(z2,Z2);
  
  ASSERT_EQ(V.size(),V2.size());
  ASSERT_EQ(Z.size(),Z2.size());
  ASSERT_EQ(V,V2);
  ASSERT_EQ(Z,Z2);
    
}

BOOST_AUTO_TEST_CASE(testCtrlObjValue){

  MtlOptEnergyTestDataModel data;
  CtrlForceEnergyTestAD ctrlFAD;
  data.initSimulator(ctrlFAD.getSimulator());
  ctrlFAD.setPenaltyCon(data.penaltyCon);
  ctrlFAD.setPartialCon(data.conF,data.conN,data.uc);
  ctrlFAD.setTotalFrames(data.T);
  ctrlFAD.setRedWarper(data.warper);
  ctrlFAD.precompute();

  ASSERT_GT (data.w.size(),0);
  const double obj_1 = data.ctrlF.fun(&data.w[0]);
  const double obj_2 = ctrlFAD.fun(&data.w[0]);
  ASSERT_EQ_TOL(obj_1,obj_2,1e-12);

}

BOOST_AUTO_TEST_CASE(testCtrlGrad){

  MtlOptEnergyTestDataModel data;
  CtrlForceEnergyTestAD ctrlFAD;
  data.initSimulator(ctrlFAD.getSimulator());
  ctrlFAD.setPenaltyCon(data.penaltyCon);
  ctrlFAD.setPartialCon(data.conF,data.conN,data.uc);
  ctrlFAD.setTotalFrames(data.T);
  ctrlFAD.setRedWarper(data.warper);
  ctrlFAD.precompute();

  ASSERT_GT (data.w.size(),0);

  VectorXd g1(ctrlFAD.dim());
  VectorXd g2(ctrlFAD.dim());

  data.ctrlF.grad(&data.w[0],&g1[0]);
  ctrlFAD.grad(&data.w[0],&g2[0]);

  ASSERT_EQ_SMALL_VEC_TOL(g1,g2,g1.size(),1e-12);

}

BOOST_AUTO_TEST_CASE(testCtrlOpt){
 
  CtrlForceEnergy ctrlF;
  
}

BOOST_AUTO_TEST_SUITE_END()
