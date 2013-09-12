#ifndef _MAMTLENERGYAD_H_
#define _MAMTLENERGYAD_H_

#include <CASADITools.h>
#include <RedRSWarperAD.h>
#include <MASimulatorAD.h>
using CasADi::SXFunction;
using namespace CASADI;
using namespace LSW_WARPING;

namespace LSW_ANI_EDITOR{

  class MAEditEnergyADItem{
	
  public:
	virtual void addEnergy(SXMatrix &energy) = 0;
  };
  typedef boost::shared_ptr<MAEditEnergyADItem> pMAEditEnergyADItem;

  class QuadricMAEditEnergyAD:public MAEditEnergyADItem{
	
  public:
	QuadricMAEditEnergyAD(const vector<SX> &variable, const VectorXd &oldValue, 
						  const double penalty=1.0f,const bool scaleByOldValue = false){

	  assert_ge(penalty,0);
	  _variable = variable;
	  _oldValue = oldValue;
	  _penalty = penalty;
	  _scaleByOldValue = scaleByOldValue;
	}
	QuadricMAEditEnergyAD(const vector<SX> &variable, const double penalty=1.0f){

	  assert_ge(penalty,0);
	  _variable = variable;
	  _oldValue = VectorXd::Zero(variable.size());
	  _penalty = penalty;
	  _scaleByOldValue = false;
	}
	QuadricMAEditEnergyAD(const SX &variable,const double oldValue, 
						  const double penalty=1.0f,const bool scaleByOldValue = false){

	  assert_ge(penalty,0);
	  _variable.resize(1);
	  _variable[0] = variable;
	  _oldValue.resize(1);
	  _oldValue[0] = oldValue;
	  _penalty = penalty;
	  _scaleByOldValue = scaleByOldValue;
	}
	virtual void addEnergy(SXMatrix &energy){

	  if(_penalty <= 0 || _variable.size() <= 0){
		return;
	  }
	  const double scale = 1.0f/(_variable.size()*_variable.size());
	  assert_eq(_variable.size(),_oldValue.size());
	  for (size_t i = 0; i < _variable.size(); ++i){
		if(_oldValue[i] != 0){
		  if(_scaleByOldValue){
			energy += (_variable[i]/_oldValue[i]-1.0)*(_variable[i]/_oldValue[i]-1.0)*_penalty*scale;
		  }else{
			energy += (_variable[i]-_oldValue[i])*(_variable[i]-_oldValue[i])*_penalty*scale;
		  }
		}else{
		  energy += _variable[i]*_variable[i]*_penalty*scale;
		}
	  }
	}

  protected:
	double _penalty;
	vector<SX> _variable;
	VectorXd _oldValue;
	bool _scaleByOldValue;
  };

  class MALinearConEnergyAD:public MAEditEnergyADItem{

  public:
	MALinearConEnergyAD(const double con_penalty,const vector<int>&g,
						const VectorXd&uc,const MatrixXd &W,
						const vector<SX>&z,const VectorXd &input_u){

	  assert_gt(con_penalty,0);
	  assert_eq((int)uc.size(),g.size()*3);
	  assert_eq(input_u.size(),W.rows());
	  
	  const int r = (int)z.size();
	  MatrixXd C(g.size()*3, r);
	  VectorXd input_uc(uc.size());
	  for (size_t k = 0; k < g.size(); ++k){
		C.block(k*3,0,3,r) = W.block(g[k]*3,0,3,r);
		input_uc.segment(k*3,3) = input_u.segment(g[k]*3,3);
	  }
	  _C = CASADI::convert(C);
	  _uc = CASADI::convert((VectorXd)uc);
	  _input_uc = CASADI::convert((VectorXd)input_uc);
	  _z.resize(r,1);
	  for (size_t i = 0; i < r; ++i){
		_z(i,0) = z[i];
	  }
	  _conPenalty = con_penalty;

	  assert_eq(_C.size1(),_uc.size1());
	  assert_eq(_z.size2(),1);
	  assert_eq(_uc.size2(),1);
	}
	virtual void addEnergy(SXMatrix &energy){

	  if(_uc.size() <=0 || _conPenalty <= 0){
		return;
	  }

	  assert_eq(_C.size2(),_z.size1());
	  const SXMatrix diff = _C.mul(_z) + _input_uc - _uc;
	  SXMatrix en = 0.0f;
	  for (int i = 0; i < diff.size1(); ++i){
		en = en + diff(i,0)*diff(i,0);
	  }
	  const double scale = 1.0f/(_uc.size()*_uc.size());
	  energy = energy + en*_conPenalty*scale;
	}

  protected:
	double _conPenalty;
	SXMatrix _z;
	SXMatrix _uc; // con disp for one frame.
	SXMatrix _C;  // con matrix of z for one frame.
	SXMatrix _input_uc;
  };

  class MAWarpConEnergyAD:public MALinearConEnergyAD{

  public:
	MAWarpConEnergyAD(const double con_penalty,const vector<int>&g,
					  const VectorXd&uc,const MatrixXd &U,
					  const vector<SX>&z,const VectorXd &input_u, 
					  pRedRSWarperAD redrs):
	  MALinearConEnergyAD(con_penalty,g,uc,U,z,input_u){
	  assert(redrs);
	  _redrs = redrs;
	  const SparseMatrix<double> &G = redrs->get_G();
	  RSCoordComp::constructWithoutWarp(G,input_u,_input_y);
	  _input_y = _redrs->getCubatureY(_input_y);
	}
	void addEnergy(SXMatrix &energy){

	  if(_uc.size() <=0 || _conPenalty <= 0){
		return;
	  }

	  assert(_redrs);
	  SXMatrix q;
	  _redrs->computeQ(_input_y,_z,q);
	  assert_eq(q.size1(),_C.size2());
	  assert_eq(q.size2(),1);
	  const SXMatrix diff = _C.mul(q) - _uc;
	  SXMatrix en = 0.0f;
	  for (int i = 0; i < diff.size1(); ++i){
		en = en + diff(i,0)*diff(i,0);
	  }
	  const double scale = 1.0f/(_uc.size()*_uc.size());
	  energy = energy + en*_conPenalty*scale;
	}

  private:
	pRedRSWarperAD _redrs;
	VectorXd _input_y;
  };

  class MAEditEnergyAD{
	
  public:
	void pushVariable(const SX &x){
	  _x.push_back(x);
	}
	void pushVariables(const vector<SX> &x){
	  _x.insert(_x.end(),x.begin(),x.end());
	}
	void pushEnergyItem(pMAEditEnergyADItem energy){
	  _energyItems.push_back(energy);
	}
	void initialize(){
	  
	  _totalEnergy = 0.0f;
	  for (size_t i = 0; i < _energyItems.size(); ++i){
		_energyItems[i]->addEnergy(_totalEnergy);
	  }
	  _totalEnergy = _totalEnergy * 0.5f;
	  _energyFun = SXFunction(_x, _totalEnergy);
	  _energyFun.init();
	}
	void clear(){
	  _x.clear();
	  _totalEnergy.clear();
	  _energyItems.clear();
	}

	const SXFunction &getEnergyFun()const{
	  return _energyFun;
	}
	SXFunction &getEnergyFun(){
	  return _energyFun;
	}
	int numVariables()const{
	  return (int)_x.size();
	}

  protected:
	SXFunction _energyFun;
	vector<SX> _x;
	SXMatrix _totalEnergy;
	vector<pMAEditEnergyADItem> _energyItems;
  };

  enum ENERGY_ITEM{CTRL_FORCE = 1, 
				   TIME_STEP = 2, 
				   FREQUENCE = 4, 
				   STIFF_DAMPING = 8, 
				   MASS_DAMPING = 16, 
				   LINEAR_CON = 32, 
				   WARP_CON = 64,
  				   STIFF_SCALE = 128};
  class AllMAEditEnergyAD{
	
  public:
	AllMAEditEnergyAD(const int energyItem=CTRL_FORCE|WARP_CON):_energyItem(energyItem){
	  
	  _penaltyEw = 1.0f; // forces.
	  _penaltyEc = 1.0f; // constraint.
	  _penaltyEh = 1.0f; // time step.
	  _penaltyEs = 1.0f; // stiffness scale.
	  _penaltyEa = 1.0f; // alpha: damping.
	  _penaltyEl = 1.0f; // lambda: eigenvalues.
	}

	void setPenaltyEw(const double penaltyEw){
	  _penaltyEw = penaltyEw;
	}
	void setPenaltyEc(const double penaltyEc){
	  _penaltyEc = penaltyEc;
	}
	void setPenaltyEh(const double penaltyEh){
	  _penaltyEh = penaltyEh;
	}
	void setPenaltyEs(const double penaltyEs){
	  _penaltyEs = penaltyEs;
	}
	void setPenaltyEa(const double penaltyEa){
	  _penaltyEa = penaltyEa;
	}
	void setPenaltyEl(const double penaltyEl){
	  _penaltyEl = penaltyEl;
	}

	void setEigenValues(const VectorXd &eigenValues){

	  _eigenValues.resize(eigenValues.size());
	  for (int i = 0; i < eigenValues.size(); ++i){
		_eigenValues[i] = eigenValues[i] > 0 ? eigenValues[i]:0; //@todo
	  }
	  _numericalSimul.setEigenValues(_eigenValues);
	  if(_energyItem & FREQUENCE){
		_eigenValues = CASADI::makeSymbolic(eigenValues.size(),"lambda");
		_frequence=pMAEditEnergyADItem(new QuadricMAEditEnergyAD(_eigenValues,eigenValues,_penaltyEl,true));
	  }
	  vector<SX> scaledLambda = _eigenValues;
	  _stiffScale = SX("stiff_scale");
	  if(_energyItem & STIFF_SCALE){
		for (size_t i = 0; i < scaledLambda.size(); ++i){
		  scaledLambda[i] *= _stiffScale;
		}
	  }
	  _symbolSimul.setEigenValues(scaledLambda);
	}
	void setTotalFrame(const int T){
	  assert_ge(T,3);
	  _T = T;
	}
	void setTimeStep(const double h){
	  assert_gt(h,0.0f);
	  _h = h;
	  _numericalSimul.setTimeStep(_h);
	  if(_energyItem & TIME_STEP){
		_h = SX("h");
		_time = pMAEditEnergyADItem(new QuadricMAEditEnergyAD(_h,h,_penaltyEh,true));
	  }
	  _symbolSimul.setTimeStep(_h);
	}
	void setStiffnessDamping(const double s){
	  VectorXd vs(reducedDim());
	  for (int i = 0; i < reducedDim(); ++i){
		vs[i] = s;
	  }
	  setStiffnessDamping(vs);
	}
 	void setMassDamping(const double s){
	  VectorXd vs(reducedDim());
	  for (int i = 0; i < reducedDim(); ++i){
		vs[i] = s;
	  }
	  setMassDamping(vs);
	}
	void setStiffnessDamping(const VectorXd &s){

	  _alphaK.resize(reducedDim());
	  for (size_t i = 0; i < _alphaK.size(); ++i){
		assert_ge(s[i],0.0f);
		_alphaK[i] = s[i];
	  }
	  _numericalSimul.setStiffnessDamping(_alphaK);
	  if(_energyItem & STIFF_DAMPING){
		assert_gt(reducedDim(),0);
		_alphaK = CASADI::makeSymbolic(reducedDim(),"alpha_k");
		_stiffDamp = pMAEditEnergyADItem(new QuadricMAEditEnergyAD(_alphaK,s,_penaltyEa,true));
	  }
	  _symbolSimul.setStiffnessDamping(_alphaK);
	}
 	void setMassDamping(const VectorXd &s){

	  assert_eq(s.size(),reducedDim());
	  _alphaM.resize(reducedDim());
	  for (size_t i = 0; i < _alphaM.size(); ++i){
		assert_ge(s[i],0.0f);
		_alphaM[i] = s[i];
	  }
	  _numericalSimul.setMassDamping(_alphaM);
	  if(_energyItem & MASS_DAMPING){
		assert_gt(reducedDim(),0);
		_alphaM = CASADI::makeSymbolic(reducedDim(),"alpha_m");
		_massDamp = pMAEditEnergyADItem(new QuadricMAEditEnergyAD(_alphaM,s,_penaltyEa,true));
	  }
	  _symbolSimul.setMassDamping(_alphaM);
	}
	void setEigenVectors(const MatrixXd &W){
	  _eigenVectors = W;
	}
	void setIntialStatus(const VectorXd &v0, const VectorXd &z0){
	  _v0 = v0;
	  _z0 = z0;
	}
	void setConstraints(const int f,const vector<int>&g,const VectorXd&uc){
	  assert_in(f,1,_T);
	  _conFrameId.push_back(f);
	  _conNodes.push_back(g);
	  _uc.push_back(uc);
	}
	void setInputU(const vector<VectorXd> &input_u){
	  _input_u = input_u;
	}
	void setRedRS(pRedRSWarperAD redrs){
	  _redrs = redrs;
	}

	void genEnergyFun(){

	  x0.resize(0);
	  xmin.resize(0);
	  
	  // forward simulate
	  const int Tw = _T-1;
	  const int r = reducedDim();
	  const vector<SX> w = CASADI::makeSymbolic(Tw*r,"w");
	  const vector<vector<SX> > z = forward(w,_symbolSimul);

	  // ctrl forces energy
	  if(_energyItem&CTRL_FORCE){
		_energy.pushVariables(w);
		_energy.pushEnergyItem(pMAEditEnergyADItem(new QuadricMAEditEnergyAD(w,_penaltyEw)));
		x0.resize(w.size(), 0);
		xmin.resize(w.size(), -numeric_limits<double>::infinity());
	  }

	  // constraint energy
	  for (size_t i = 0; i < _conFrameId.size(); ++i){

		assert_in(_conFrameId[i],0,z.size());
		assert_eq(_input_u.size(),z.size()+1);
		const vector<SX> &zi = z[_conFrameId[i]-1];
		const VectorXd &input_u = _input_u[_conFrameId[i]];
		assert_eq((int)zi.size(),reducedDim());
		pMAEditEnergyADItem con;
		if(_energyItem & LINEAR_CON){
		  con = pMAEditEnergyADItem(new MALinearConEnergyAD(_penaltyEc,_conNodes[i],_uc[i],_eigenVectors,zi,input_u));
		}else if(_energyItem & WARP_CON){
		  const MatrixXd &U = _redrs->getBasisMat();
		  con = pMAEditEnergyADItem(new MAWarpConEnergyAD(_penaltyEc,_conNodes[i],_uc[i],U,zi,input_u,_redrs));
		}
		if(con){
		  _energy.pushEnergyItem(con);
		}
	  }

	  // material energies
	  if(_energyItem & TIME_STEP){
		_energy.pushVariable(_h);
		_energy.pushEnergyItem(pMAEditEnergyADItem(_time));
		x0.resize(x0.size()+1, _numericalSimul.getTimeStep().getValue());
		xmin.resize(xmin.size()+1, 0);
	  }
	  if(_energyItem & STIFF_SCALE){
		_energy.pushVariable(_stiffScale);
		_energy.pushEnergyItem(pMAEditEnergyADItem(new QuadricMAEditEnergyAD(_stiffScale,1.0f,_penaltyEs,false)));
		x0.resize(x0.size()+1, 1.0);
		xmin.resize(xmin.size()+1, 0);
	  }
	  if(_energyItem & STIFF_DAMPING){
		_energy.pushVariables(_alphaK);
		_energy.pushEnergyItem(_stiffDamp);
		x0.resize(x0.size()+_alphaK.size(), _numericalSimul.getStiffnessDamping()[0].getValue());
		xmin.resize(xmin.size()+_alphaK.size(), 0);
	  }
	  if(_energyItem & MASS_DAMPING){
		_energy.pushVariables(_alphaM);
		_energy.pushEnergyItem(_massDamp);
		x0.resize(x0.size()+_alphaM.size(), _numericalSimul.getMassDamping()[0].getValue());
		xmin.resize(xmin.size()+_alphaM.size(), 0);
	  }
	  if(_energyItem & FREQUENCE){

		_energy.pushVariables(_eigenValues);
		_energy.pushEnergyItem(_frequence);
		for (size_t i = 0; i < _eigenValues.size(); ++i){
		  x0.push_back(_numericalSimul.getEigenValues()[i].getValue());
		}
		xmin.resize(xmin.size()+_eigenValues.size(),0);
	  }

	  // initialize 
	  _energy.initialize();
	}
	const vector<double> &getInitX()const{
	  return x0;
	}
	const vector<double> &getBoundX()const{
	  return xmin;
	}

	const SXFunction &getEnergyFun()const{
	  return _energy.getEnergyFun();
	}
	SXFunction &getEnergyFun(){
	  return _energy.getEnergyFun();
	}
	int numVariables()const{
	  return _energy.numVariables();
	}
	int reducedDim()const{
	  return _eigenValues.size();
	}
	int totalFrames()const{
	  return _T;
	}
	vector<VectorXd> getResultZ(const VectorXd &w){

	  vector<SX> ws(w.size());
	  for (int i = 0; i < w.size(); ++i){
		ws[i] = w[i];
	  }
	  const vector<vector<SX> > zs = forward(ws,_numericalSimul);
	  vector<VectorXd> z;
	  z.push_back(_z0);
	  for (size_t f = 0; f < zs.size(); ++f){
		VectorXd zf(zs[f].size());
		for (int i = 0; i < zs[f].size(); ++i){
		  zf[i] = zs[f][i].getValue();
		}
		z.push_back(zf);
	  }
	  return z;
	}
	int getEnergyItem()const{
	  return _energyItem;
	}

  protected:
	const vector<vector<SX> > forward(const vector<SX> &w,MASimulatorAD &simulator)const{

	  vector<vector<SX> > z;
	  const int Tw = _T-1;
	  const int r = reducedDim();
	  assert_eq((int)w.size(),Tw*r);

	  simulator.setIntialStatus(_v0,_z0);
	  for (int i = 0; i < Tw; ++i){
		const vector<SX> wi(w.begin()+i*r,w.begin()+(i+1)*r);
		simulator.forward(wi);
		z.push_back(simulator.getZ());
	  }
	  return z;
	}

  private:
	const int _energyItem;
	MASimulatorAD _symbolSimul;
	MASimulatorAD _numericalSimul;
	MAEditEnergyAD _energy;
	pMAEditEnergyADItem _time;
	pMAEditEnergyADItem _stiffDamp;
	pMAEditEnergyADItem _massDamp;
	pMAEditEnergyADItem _frequence;
	pRedRSWarperAD _redrs;

	double _penaltyEw; // forces.
	double _penaltyEc; // constraint.
	double _penaltyEh; // time step.
	double _penaltyEs; // stiffness scale.
	double _penaltyEa; // alpha: damping.
	double _penaltyEl; // lambda: eigenvalues.

	int _T; 
	SX _h; 
	SX _stiffScale;
	vector<SX> _alphaK;
	vector<SX> _alphaM;
	vector<SX> _eigenValues;
	MatrixXd _eigenVectors;
	VectorXd _v0;
	VectorXd _z0;
	vector<VectorXd> _input_u;

	vector<int> _conFrameId;
	vector<vector<int> > _conNodes;
	vector<VectorXd> _uc;

	vector<double> x0; // initial value for var.
	vector<double> xmin; // boundary for var.
  };
}

#endif /* _MAMTLENERGYAD_H_ */
