#ifndef _FULLCONMAT_H_
#define _FULLCONMAT_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Sparse>
#include <assertext.h>
using namespace Eigen;
using namespace std;

namespace LSW_WARPING{
  
  /**
   * @class FullConMat Compute the constraint matrix when the index of the
   * constrained nodes provided, which are 0-indexed. The constraint matrix is
   * used to fixed all of the provided nodes fully, not just the barycenter.
   * 
   */
  class FullConMat{
	
  public:
	static void compute(const vector<int> &con_nodes, SparseMatrix<double> &C, const int total_node_num){
  
	  assert_ge(total_node_num,0);
	  assert_ge(con_nodes.size(),0);

	  C.resize(con_nodes.size()*3, total_node_num*3);
	  C.reserve(con_nodes.size()*3);
	  for (int i = 0; i < (int)con_nodes.size(); ++i){
		insert_I3x3(C, con_nodes[i], i);
	  }
	}
	static void insert_I3x3(SparseMatrix<double> &m, const int node_id, const int c_index){
  
	  assert_ge(node_id,0);
	  assert_ge(m.cols(),3);
	  assert_ge(m.rows(),3);
	  assert_gt(m.cols(),node_id*3+2);
	  assert_gt(m.rows(),c_index*3+2);

	  m.insert(c_index*3+0,node_id*3+0) = 1;
	  m.insert(c_index*3+1,node_id*3+1) = 1;
	  m.insert(c_index*3+2,node_id*3+2) = 1;
	}
	
  };
  
  typedef boost::shared_ptr<FullConMat> pFullConMat;
  
}//end of namespace

#endif /*_FULLCONMAT_H_*/
