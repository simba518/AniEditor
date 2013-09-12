#ifndef _REDBARYCENCONM_H_
#define _REDBARYCENCONM_H_

#include <assert.h>
#include <vector>
#include <set>
using namespace std;

#include <eigen3/Eigen/Dense>
using namespace Eigen;

namespace LSW_SIM{

  /**
   * @class RedBaryCenConM
   * @brief compute the constraint matrixes for barycenter constraint of 
   * reduced stvk simulation method.
   */
  class RedBaryCenConM{
	
  public:
	static const MatrixXd &computeConM(const set<int> &con_nodes, 
									   const MatrixXd &U,
									   MatrixXd &con_m);

	void computeConM(const vector<set<int> >&con_nodes,const MatrixXd &U);

	const MatrixXd &getConM(int i)const{
	  assert(i>=0 && i<(int)con_Ms.size());
	  return con_Ms[i];
	}

	const vector<MatrixXd > &getConM()const{
	  return con_Ms;
	}

  protected:
	
  private:
	vector<MatrixXd > con_Ms;

  };
  
}//end of namespace

#endif /*_REDBARYCENCONM_H_*/
