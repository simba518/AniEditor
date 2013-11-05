#ifndef _BASEOPTFUN_H_
#define _BASEOPTFUN_H_

#include <boost/shared_ptr.hpp>

namespace LSW_ANI_EDITOR{

  class BaseOptFun{
	
  public:
	virtual int dim()const = 0;
	virtual void init(double *x,const int n) = 0;
	virtual void bounds(double *x_l,double *x_u,const int n){
	  for (int i = 0; i < n; ++i){
		x_l[i] = -2e19;
		x_u[i] = 2e19;
	  }
	}
	virtual double fun(const double *x) = 0;
	virtual void grad(const double *x,double *g) = 0;
	virtual void setRlst(const double *x, const double objValue) = 0;
  };
  typedef boost::shared_ptr<BaseOptFun> pBaseOptFun;
  
}//end of namespace

#endif /*_BASEOPTFUN_H_*/
