#ifndef _SPACETIMEHESSIANAUTODIFF_H_
#define _SPACETIMEHESSIANAUTODIFF_H_

#include <vector>
using namespace std;

#include <boost/shared_ptr.hpp>
#include <BaseSpaceTimeHessian.h>
#include <symbolic/sx/sx.hpp>
#include <symbolic/fx/sx_function.hpp>

namespace IEDS{
  
  /**
   * @class SpaceTimeHessianAutoDiff compute the hessian of the spacetime
   * function using auto difference approach.
   * 
   * @see BaseSpaceTimeHessian.h
   */
  class SpaceTimeHessianAutoDiff: public BaseSpaceTimeHessian{
	
  public:
	SpaceTimeHessianAutoDiff();
	bool compute_HTriplet(vector<E_Triplet> &H_triplet);

  protected:
	typedef vector<CasADi::SX>::iterator vsx_it;
	void make_symbolic(vsx_it begin, vsx_it end, const string name){

	  int counter = 0;
	  for (vsx_it i = begin; i != end; ++i){
		(*i) = CasADi::SX(name+string("_")+boost::lexical_cast<string>(counter));
		counter ++;
	  }
	}
	void makeSymbol_H();
	void calculate_H(vector<E_Triplet> &H_triplet);

  private:
	CasADi::SXMatrix H_SX;

  };
  
  typedef boost::shared_ptr<SpaceTimeHessianAutoDiff> pSpaceTimeHessianAutoDiff;
  
}//end of namespace

#endif /*_SPACETIMEHESSIANAUTODIFF_H_*/
