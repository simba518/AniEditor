#ifndef _BARYCENTERDRAGGERVEC_H_
#define _BARYCENTERDRAGGERVEC_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
using Eigen::VectorXd;

#include <BarycenterDragger.h>

namespace LSW_SIM{
  
  /**
   * @class BarycenterDraggerVec a class that warps the class BarycenterDragger.
   * In particular, use VectorXd instead of "double *" in the setting up
   * functions for safety.
   */
  class BarycenterDraggerVec: public BarycenterDragger{
	
  public:
	void setConGroups(const vector<set<int> > &groups, const VectorXd &vol_u);
	void stopDrag();
	VectorXd getUcVec()const;

  protected:
	bool checkValidDation(const vector<set<int> > &g, const VectorXd &u)const;
	
  private:
	
  };
  
  typedef boost::shared_ptr<BarycenterDraggerVec> pBarycenterDraggerVec;
  
}//end of namespace

#endif /*BarycenterDraggerVec*/
