#include <boost/foreach.hpp>

#include "FullBaryCenConM.h"
using namespace LSW_SIM;

#include "RedCenConM.h"
using namespace LSW_SIM;

MatrixXd &RedCenConM::computeConM(const vector<set<int> > &con_nodes,
								  const MatrixXd &U, MatrixXd &con_m){

  const int total_node_num = U.rows()/3;
  const int r = U.cols();
  con_m.resize(3*con_nodes.size(),r);

  int i = 0;
  BOOST_FOREACH(const set<int> &group, con_nodes){

	SparseMatrixXd sparse_con_m;
	FullBaryCenConM::computeConM(group,total_node_num,sparse_con_m);
	con_m.block(i*3,0,3,r) = sparse_con_m*U;
	++i;
  }
  return con_m;
}
