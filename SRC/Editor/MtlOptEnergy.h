#ifndef _MTLOPTENERGY_H_
#define _MTLOPTENERGY_H_

#include <BaseMtlOptEnergy.h>

namespace LSW_ANI_EDITOR{
  
  class CtrlForceEnergy:public BaseCtrlForceEnergy,public BaseOptFun{
	
  public:
	int dim()const{
	  
	}
	void init(double *x,const int n){
	  assert_gt(n,0);
	  assert_eq(n,dim());
	}
	double fun(const double *x){}
	void grad(const double *x,double *g){}
	void setRlst(const double *x, const double objValue){}

  protected:
	void precompute(){
	  
	}
	
  private:
	
  };
  typedef boost::shared_ptr<CtrlForceEnergy> pCtrlForceEnergy;

  class MtlOptEnergy: public BaseMtlOptEnergy{
	
  public:
	double fun(const double *x){
	  computeKD(x,_K,_D);
	}
	void grad(const double *x,double *g){}

  protected:
	
  };
  typedef boost::shared_ptr<MtlOptEnergy> pMtlOptEnergy;

}//end of namespace

#endif /* _MTLOPTENERGY_H_ */
