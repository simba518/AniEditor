#ifndef _UNWARPERAD_H_
#define _UNWARPERAD_H_

#include <WarperAD.h>

namespace LSW_WARPING{
  
  class PosWarpEnergyAD{
  
  public:
	static SXMatrix compute(const vector<SXMatrix> &u, 
							const vector<SXMatrix> &w, 
							const vector<SXMatrix> &p){

	  assert_eq (u.size(),w.size());
	  assert_eq (u.size(),p.size());
	  SXMatrix energy = 0.0f;
	  for (size_t i = 0; i < w.size(); ++i){
		const SXMatrix _R = WarperAD::warpMat(w[i]);
		const SXMatrix warpedU = mul(_R,p[i]);
		energy += inner_prod(warpedU-u[i],warpedU-u[i]);
	  }
	  return energy;
	}
  };

  class GradWarpEnergyAD{
	
  public:
	static SXMatrix compute(pTetMesh_const rest, const vector<SXMatrix> &u,
							const vector<SXMatrix> &w, const vector<SXMatrix> &p){

	  SparseMatrix<double> G;
	  DefGradOperator::compute(rest,G);
	  SXMatrix sG;
	  CASADI::convert(G,sG);

	  SXMatrix u_vec(u.size()*3,1);
	  for (size_t i = 0; i < u.size(); ++i){
		assert_eq(u[i].size1(),3);
		assert_eq(u[i].size2(),1);
		u_vec(i*3+0) = u[i](0);
		u_vec(i*3+1) = u[i](1);
		u_vec(i*3+2) = u[i](2);
	  }
	  const SXMatrix gu = mul(sG,u_vec);

	  SXMatrix warpedU(u.size()*3,1);
	  SXMatrix u3(3,1);
	  for (size_t i = 0; i < w.size(); ++i){
		
		const SXMatrix _R = WarperAD::warpMat(w[i]);
		const SXMatrix warpedU3 = mul(_R,p[i]);
		warpedU(i*3+0) = warpedU3(0);
		warpedU(i*3+1) = warpedU3(1);
		warpedU(i*3+2) = warpedU3(2);
	  }

	  const SXMatrix gru = mul(sG,warpedU);

	  /// @todo add volume as weights for each element
	  return inner_prod(gru-gu,gru-gu);
	}

  };

  class RotWarpEnergyAD{
	
  public:
	static SXMatrix compute(pTetMesh_const rest, 
							const vector<SXMatrix> &w, 
							const vector<SXMatrix> &p){

	  SXMatrix pv(p.size()*3,1);
	  for (size_t i = 0; i < p.size(); ++i){
	  	assert_eq(p[i].size(),3);
	  	pv(i*3+0) = p[i](0);
	  	pv(i*3+1) = p[i](1);
	  	pv(i*3+2) = p[i](2);
	  }
	  
	  vector<SXMatrix> nodeW;
	  WarperAD::computeRotVector(rest,pv,nodeW);
	  
	  // compute energy
	  SXMatrix energy = 0.0f;
	  for (size_t i = 0; i < w.size(); ++i){
		energy += inner_prod(w[i]-nodeW[i], w[i]-nodeW[i]);
	  }
	  return energy;
	}
	
  };

  class UnwarpEnergyAD{
	
  public:
	UnwarpEnergyAD(pTetMesh_const rest, double a, double b, double g):
	  _restMesh(rest),_alpha(a),_beta(b),_gamma(g){

	  assert(rest != NULL);
	  assert_ge(_alpha,0.0f);
	  assert_ge(_beta,0.0f);
	  assert_ge(_gamma,0.0f);
	  assert_in(_gamma+_alpha+_beta,1.0f-1e-3,1.0f+1e-3);
	  precompute();
	}
	// compute _energyFun_wp
	void genEnergyFun(const VectorXd &u){
	  
	  assert_gt(u.size(),0);
	  assert_eq(u.size()*2, getDim());
	  const vector<SX> w = CASADI::makeSymbolic(getDim()/2,"w");
	  const vector<SX> p = CASADI::makeSymbolic(getDim()/2,"p");
	  const vector<SX> w_p = CASADI::connect(w,p);
	  vector<SX> us(u.size());
	  for (size_t i = 0; i < u.size(); ++i){
		us[i] = u[i];
	  }
	  const vector<SX> u_w_p = CASADI::connect(us,w_p);
	  _energyFun_wp =SXFunction(w_p,_energyFun_uwp.eval(u_w_p));
	}
	const SXFunction &getEnergyFun()const{
	  return _energyFun_wp;
	}
	int getDim()const{
	  if (_restMesh)
		return _restMesh->nodes().size()*6;
	  return 0;
	}

  protected:
	// intialize _energyFun_uwp
	void precompute(){
	  
	  const vector<SXMatrix> u = CASADI::Sx2Mat(CASADI::makeSymbolic(getDim()/2,"u"),getDim()/6);
	  const vector<SXMatrix> w = CASADI::Sx2Mat(CASADI::makeSymbolic(getDim()/2,"w"),getDim()/6);
	  const vector<SXMatrix> p = CASADI::Sx2Mat(CASADI::makeSymbolic(getDim()/2,"p"),getDim()/6);
	  const SXMatrix energy 
	  	= _alpha*PosWarpEnergyAD::compute(u,w,p) +
	  	_beta*GradWarpEnergyAD::compute(_restMesh,u,w,p) +
	  	_gamma*RotWarpEnergyAD::compute(_restMesh,w,p);
	  const vector<SX> w_p = CASADI::connect( CASADI::Mat2Sx(w),CASADI::Mat2Sx(p) );
	  const vector<SX> u_w_p = CASADI::connect( CASADI::Mat2Sx(u),w_p);
	  _energyFun_uwp = SXFunction(u_w_p,energy);
	  _energyFun_uwp.init();
	}
	
  private:
	SXFunction _energyFun_uwp; // f(u,w,p)
	SXFunction _energyFun_wp;  // f(w,p)
	pTetMesh_const _restMesh;
	double _alpha;
	double _beta;
	double _gamma;
  };

  typedef boost::shared_ptr<UnwarpEnergyAD> pUnwarpEnergyAD;

}
#endif /* _UNWARPERAD_H_ */
