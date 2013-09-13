#ifndef _FULLBARYCENCONM_H_
#define _FULLBARYCENCONM_H_

#include <assert.h>
#include <vector>
#include <set>
using namespace std;

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <eigen3/Eigen/Sparse>
using namespace Eigen;

typedef SparseMatrix<double> SparseMatrixXd;

namespace LSW_SIM{

  typedef vector<Eigen::Triplet<double> > VecT;
  
  class FullBaryCenConM{
	
  public:
	static const SparseMatrixXd &computeConM(const set<int> &con_nodes, 
											 int total_node_num,
											 SparseMatrixXd &con_m);

	static const SparseMatrixXd &computeConM(const vector<set<int> > &con_nodes, 
											 int total_node_num,
											 SparseMatrixXd &con_m);

	static const VecT &computeConM(const vector<set<int> > &con_nodes, 
								   int total_node_num,
								   VecT &con_M_triplet);

	void computeConM(const vector<set<int> >&con_nodes, 
					 int total_node_num);

	const SparseMatrixXd &getConM(int i)const{
	  assert(i>=0 && i<(int)con_Ms.size());
	  return con_Ms[i];
	}

  protected:
	

  private:
	vector<SparseMatrixXd > con_Ms;
  };
  
}//end of namespace

#endif /*_FULLBARYCENCONM_H_*/
