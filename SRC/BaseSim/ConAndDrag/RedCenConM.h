
#ifndef _REDCENCONM_H_
#define _REDCENCONM_H_

#include <set>
#include <vector>
using namespace std;

#include <boost/shared_ptr.hpp>

#include <eigen3/Eigen/Dense>
using Eigen::MatrixXd;

namespace LSW_SIM{
  
  /**
   * @class RedCenConM compute the constraint matrix C in subspace for Cu = uc.
   * When there are 3 constrained groups, we have
   *     |C1|           |uc1|
   * C = |C2|  and uc = |uc2|
   *     |C3|           |uc3|
   * 
   * Here, uci is the displacements of barycenter of i-th constrained groups,
   * and Ci is the corresponding constraint matrix. uci is a 3x1 vector, and Ci
   * is 3x(node_number*3).
   * 
   */
  class RedCenConM{
	
  public:
	/**
	 * con_nodes is the nodes' ids of the constrained groups ,U is the nonlinear
	 * basis, and con_m is the result constraint matrix.
	 */
	static MatrixXd &computeConM(const vector<set<int> > &con_nodes,
								 const MatrixXd &U, MatrixXd &con_m);
	
  };
  
  typedef boost::shared_ptr<RedCenConM> pRedCenConM;
  
}//end of namespace

#endif /*_REDCENCONM_H_*/
