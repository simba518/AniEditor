#include "MtlOpt.h"

void RedSpaceTimeEnergyAD::insertKeyframes(MatrixXd &Z)const{

  const int r = reducedDim();
  const int T = Z.cols()+_keyZ.cols();
  assert_eq(Z.rows(),r);
  assert_eq(T,_T);
  
  const MatrixXd oldZ = Z;
  Z.resize(r,T);
  Z.setZero();
  int xpos = 0;
  for (size_t i = 0; i < T; ++i){
  	const int k = isKeyframe(i);
  	if(k >= 0){
  	  Z.col(i) = _keyZ.col(k);
  	}else{
  	  Z.col(i) = oldZ.col(xpos);
  	  xpos++;
  	}
  }
}

void MtlOptimizer::optimize(){

  this->optimizationBegin();

  _energy.assembleEnergy();
  const SXMatrix &E = _energy.getEnergy();
  const VSX &x = this->getVariable();

  _fun = CasADi::SXFunction(x,E);
  _solver = CasADi::IpoptSolver(_fun);
  _solver.setOption("generate_hessian",true);
  _solver.init();

  VectorXd x0;
  this->getInitValue(x0,x.size());
  assert_gt(x0.size(),0);
  assert_eq(x0.size(),x.size());
  _solver.setInput(&x0[0],CasADi::NLP_X_INIT);

  vector<double> lower;
  if(getLowBound(lower)){
	assert_eq(lower.size(), x0.size());
	_solver.setInput(lower,CasADi::NLP_LBX);
  }

  _solver.solve();

  std::vector<double> vx(x.size());
  _solver.getOutput(vx,CasADi::NLP_X_OPT);
  _rlst.resize(vx.size());
  for (size_t i = 0; i < vx.size(); ++i)
	_rlst[i] = vx[i];

  this->optimizationEnd();
}
  
void MtlOptimizer::resetEnergy(const bool useAllZ){

  const int T = _model.T;
  _energy.setT(T);
  _energy.setTimestep(_model.h);
  _energy.setDamping(_model.D);
  _energy.setK(_model.K);
  _energy.setKeyframes(_model.keyZ,_model.keyId);
	
  if(useAllZ){
	assert_eq(_model.Z.cols(), T);
	assert_eq(_model.Z.rows(), _energy.reducedDim());
	vector<int> kfid(T);
	for (int i = 0; i < T; ++i)
	  kfid[i] = i;
	_energy.setKeyframes(_model.Z,kfid);
  }
}

bool MtlOptimizer::getInitValue(VectorXd &x0, const int size)const{

  if(_rlst.size() == size){
	x0 = _rlst;
  }else {
	x0.resize(size);
	x0.setZero();
  }
  return true;
}
  
void MtlOptimizer::produceSymetricMat(const string name,const int dim,VSX&s,SXMatrix &SM){
	
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

int MtlOptimizer::symIndex(int r,int c){
  // 0
  // 1,2
  // 3,4,5
  // 6,7,8,9
  const int rr = r>=c?r:c;
  const int cc = r>=c?c:r;
  const int n = (rr+1)*(rr+2)/2-1-(rr-cc);
  return n;
}
