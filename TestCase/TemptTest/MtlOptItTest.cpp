#include "MtlOptModel.h"

class MtlEnergyZ{
  
public:
  MtlEnergyZ(const MtlOptModel &model):_model(model){
	model.initMtlOpt(_energy);
  }
  void setMtl(const MatrixXd &D,const MatrixXd &K){
	_energy.setDamping(D);
	_energy.setK(K);
  }
  void optimize(){

	_energy.assembleEnergy();
	const SXMatrix &E = _energy.getEnergy();
	const VSX &z = _energy.getVarZ();

	CasADi::SXFunction fun(z,E);
	CasADi::IpoptSolver solver(fun);
	solver.setOption("generate_hessian",true);
	solver.init();
	if(_zInit.size() != z.size()){
	  _zInit.resize(z.size());
	  _zInit.setZero();
	}
	solver.setInput(&_zInit[0],CasADi::NLP_X_INIT);
	solver.solve();

	vector<double> vx(z.size());
	solver.getOutput(vx,CasADi::NLP_X_OPT);
	assert_eq(vx.size(),_zInit.size());
	for (size_t i = 0; i < vx.size(); ++i)
	  _zInit[i] = vx[i];

	MatrixXd kzM = _model.Kz;
	const VectorXd kz = Map<VectorXd>(&(kzM(0,0)),kzM.size());
	const int r = _model.redDim();
	MatrixXd mz = MtlOptModel::assembleFullZ(_zInit,kz,_model.Kid,r);
	_z = Map<VectorXd>(&mz(0,0),mz.size());
  }
  const VectorXd &getZ()const{
	return _z;
  }

private:
  const MtlOptModel &_model;
  RedSpaceTimeEnergyAD _energy;
  VectorXd _zInit;
  VectorXd _z;
};

class MtlEnergyAtA{
  
public:
  MtlEnergyAtA(const MtlOptModel &model):_model(model){
	model.initMtlOpt(_energy);
  }
  void setZ(VectorXd z){

	const int T = _energy.getT();
	const int r = _energy.reducedDim();
	ASSERT_EQ(T*r,z.size());
	vector<int> kfid(T);
	for (int i = 0; i < T; ++i)
	  kfid[i] = i;
	const MatrixXd Kz = Map<MatrixXd>(&z[0],r,T);
	_energy.setKeyframes(Kz,kfid);
  }
  void optimize(){

	_kdVar.clear();
	this->setK();
	this->setD();
	this->setKDBound();
	if(_KDInit.size() != _kdVar.size()){
	  this->setInitKD();
	}
	assert_eq(_KDInit.size(),_kdVar.size());

	_energy.assembleEnergy();
	const SXMatrix &E = _energy.getEnergy();
	assert_eq(_energy.getVarZ().size(),0);

	CasADi::SXFunction fun(_kdVar,E);
	CasADi::IpoptSolver solver(fun);
	solver.setOption("generate_hessian",true);
	solver.init();

	solver.setInput(&_KDInit[0],CasADi::NLP_X_INIT);
	if(_KDLowB.size() == _kdVar.size()){
	  solver.setInput(_KDLowB,CasADi::NLP_LBX);
	}

	solver.solve();
	
	vector<double> vx(_kdVar.size());
	solver.getOutput(vx,CasADi::NLP_X_OPT);
	for (size_t i = 0; i < vx.size(); ++i)
	  _KDInit[i] = vx[i];
	
	this->computeK();
	this->computeD();
  }
  const MatrixXd& getK()const{
	return _K;
  }
  const MatrixXd& getD()const{
	return _D;
  }

protected:
  virtual void setK(){

	const int r = _model.redDim();
	const VSX ax = makeSymbolic(r*r,"a");
	const SXMatrix A = convert(ax,r);
	const SXMatrix K = CasADi::trans(A).mul(A);
	_energy.setK(K);
	_kdVar.insert(_kdVar.end(),ax.begin(), ax.end());
  }
  virtual void setD(){
	const int r = _model.redDim();
	const VSX damping = makeSymbolic(r,"d");
	_energy.setDamping(damping);
	_kdVar.insert(_kdVar.end(),damping.begin(), damping.end());
  }
  virtual void setInitKD(){
	const int r = _model.redDim();
	_KDInit.resize(r*r+r);
	MatrixXd K = _model.lambda.asDiagonal();
	for (int i = 0; i < r; ++i)
	  K(i,i) = sqrt(K(i,i));
	_KDInit.segment(0,r*r) = Map<VectorXd>(&K(0,0),r*r);
	_KDInit.segment(r*r,r) = _model.alphaM*VectorXd::Ones(r)+_model.alphaK*_model.lambda;
  }
  virtual void setKDBound(){
	// const int r = _model.redDim();
	// _KDLowB = vector<double> (r*r+r,0);
  }
  virtual void computeK(){

	const int r = _energy.reducedDim();
	const MatrixXd A = Map<MatrixXd>(&_KDInit[0],r,r);
	_K = A.transpose()*A;
	cout<< endl << endl << "************** K = " << _K << endl;
	EigenSolver<MatrixXd> eigenK(_K);
	cout<< endl << endl<< "eigen(K): " << eigenK.eigenvalues().transpose() << endl << endl;
  }
  virtual void computeD(){
	
	const int r = _energy.reducedDim();
	const VectorXd d = Map<VectorXd>(&_KDInit[r*r],r);
	_D = d.asDiagonal();
  }

protected:
  const MtlOptModel &_model;
  RedSpaceTimeEnergyAD _energy;
  VectorXd _KDInit;
  vector<double> _KDLowB;
  MatrixXd _K;
  MatrixXd _D;
  VSX _kdVar;
};

class MtlEnergyAtADenseD:public MtlEnergyAtA{
  
public:
  MtlEnergyAtADenseD(const MtlOptModel &model):MtlEnergyAtA(model){};

protected:
  virtual void setD(){

	const int r = _model.redDim();
	const VSX damping = makeSymbolic(r*r,"d");
	const SXMatrix D = convert(damping,r);
	_energy.setDamping(D);
	_kdVar.insert(_kdVar.end(),damping.begin(), damping.end());
  }
  virtual void setInitKD(){

	const int r = _model.redDim();
	_KDInit.resize(r*r*2);
	MatrixXd K = _model.lambda.asDiagonal();
	for (int i = 0; i < r; ++i)
	  K(i,i) = sqrt(K(i,i));
	const double ak = _model.alphaK;
	const double am = _model.alphaM;
	MatrixXd D = MatrixXd::Identity(r,r)*am + K*ak;

	_KDInit.segment(0,r*r) = Map<VectorXd>(&K(0,0),r*r);
	_KDInit.segment(r*r,r*r) = Map<VectorXd>(&D(0,0),r*r);
  }
  virtual void setKDBound(){
	// const int r = _model.redDim();
	// _KDLowB = vector<double> (r*r+r,0);
  }
  virtual void computeD(){
	const int r = _energy.reducedDim();
	_D = Map<MatrixXd>(&_KDInit[r*r],r,r);
  }
};

class MtlEnergyK:public MtlEnergyAtA{
  
public:
  MtlEnergyK(const MtlOptModel &model):MtlEnergyAtA(model){}

  static void produceSymetricMat(const string name,const int dim,vector<CasADi::SX> &s,CasADi::SXMatrix &SM){
	
	assert_gt(dim,0);
	assert_gt(name.size(),0);
	s = makeSymbolic(dim*(1+dim)/2,name);

	SM.makeDense(dim, dim, 0.0f);
	for (int i = 0; i < SM.size1(); ++i){
	  for (int j = 0; j < SM.size2(); ++j){
		const int n = symIndex(i,j);
		assert_in(n,0,s.size()-1);
		SM.elem(i,j) = s[n];
	  }
	}
  }
  static int symIndex(int r,int c){
	// 0
	// 1,2
	// 3,4,5
	// 6,7,8,9
	const int rr = r>=c?r:c;
	const int cc = r>=c?c:r;
	const int n = (rr+1)*(rr+2)/2-1-(rr-cc);
	return n;
  }

protected:
  virtual void setK(){

	const int r = _model.redDim();
	const int n = r*(r+1)/2;
	VSX ak;
	SXMatrix K;
	produceSymetricMat("k",r,ak,K);
	assert_eq(ak.size(),n);
	_energy.setK(K);
	_kdVar.insert(_kdVar.end(),ak.begin(), ak.end());
  }
  virtual void setInitKD(){

	const int r = _model.redDim();
	const int n = r*(r+1)/2;
	_KDInit.resize(n+r);
	MatrixXd K = _model.lambda.asDiagonal();
	for (int i = 0; i < r; ++i){
	  for (int j = 0; j < r; ++j){
		const int k = symIndex(i,j);
		_KDInit[k] = K(i,j);
	  }
	}
	_KDInit.segment(n,r) = _model.alphaM*VectorXd::Ones(r)+_model.alphaK*_model.lambda;
  }
  virtual void computeK(){
	const int r = _energy.reducedDim();
	const int n = r*(r+1)/2;
	_K.resize(r,r);
	for (int i = 0; i < r; ++i){
	  for (int j = 0; j < r; ++j){
		_K(i,j) = _KDInit( symIndex(i,j) );
	  }
	}
	cout << "************** K = " << _K << endl;
  }
  virtual void computeD(){
	const int r = _energy.reducedDim();
	const int n = r*(r+1)/2;
	const VectorXd d = Map<VectorXd>(&_KDInit[n],r);
	_D = d.asDiagonal();
  }
};

BOOST_AUTO_TEST_SUITE(MtlOptTest)

BOOST_AUTO_TEST_CASE(Opt_AtA){

  MtlOptModel model;
  model.produceSimRlst();
  // model.extrangeKeyframes();

  MtlEnergyZ optZ(model);
  MtlEnergyAtA optK(model);
  
  // const int r = model.redDim();
  // optZ.setMtl(MatrixXd::Zero(r,r),MatrixXd::Zero(r,r));
  for (int i = 0; i < 20; ++i){

	optZ.optimize();
	optK.setZ(optZ.getZ());
	optK.optimize();
	optZ.setMtl(optK.getD(),optK.getK());
	cout << "ITERATION: " << i << "---------------------------------------" << endl;
  }
}

BOOST_AUTO_TEST_CASE(Opt_K){

  MtlOptModel model;
  model.produceSimRlst();
  // model.extrangeKeyframes();

  MtlEnergyZ optZ(model);
  MtlEnergyK optK(model);
  
  // const int r = model.redDim();
  // optZ.setMtl(MatrixXd::Zero(r,r),MatrixXd::Zero(r,r));
  for (int i = 0; i < 20; ++i){

	optZ.optimize();
	optK.setZ(optZ.getZ());
	optK.optimize();
	optZ.setMtl(optK.getD(),optK.getK());
	cout << "ITERATION: " << i << "---------------------------------------" << endl;
  }
}

BOOST_AUTO_TEST_CASE(Opt_AtA_denseD){

  MtlOptModel model;
  model.produceSimRlst();
  model.extrangeKeyframes();

  MtlEnergyZ optZ(model);
  MtlEnergyAtADenseD optK(model);
  
  // const int r = model.redDim();
  // optZ.setMtl(MatrixXd::Zero(r,r),MatrixXd::Zero(r,r));
  for (int i = 0; i < 30; ++i){

	optZ.optimize();
	optK.setZ(optZ.getZ());
	optK.optimize();
	optZ.setMtl(optK.getD(),optK.getK());
	cout << "ITERATION: " << i << "---------------------------------------" << endl;
  }

  VectorXd z = optZ.getZ();
  assert_eq(z.size(),model.T*model.redDim());
  const MatrixXd Z = Map<MatrixXd>(&z[0],model.redDim(),model.T);
  model.saveMesh(Z,"Opt_AtA_denseD");
  model.saveMesh(model.Z,"input");
}

BOOST_AUTO_TEST_SUITE_END()
