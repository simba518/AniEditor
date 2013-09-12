#include "RedBaryCenConM.h"
#include "FullBaryCenConM.h"

using namespace LSW_SIM;

const MatrixXd &RedBaryCenConM::computeConM(const set<int> &con_nodes,const MatrixXd &U,MatrixXd &con_m){

  SparseMatrixXd sparse_con_m;
  int total_node_num = U.rows() / 3;
  FullBaryCenConM::computeConM(con_nodes,total_node_num,sparse_con_m);
  con_m = sparse_con_m * U;
  return con_m;
}

void RedBaryCenConM::computeConM(const vector<set<int> >&con_nodes, const MatrixXd &U){

  con_Ms.resize(con_nodes.size());
  vector<set<int> >::const_iterator it = con_nodes.begin();
  for( int i = 0; it != con_nodes.end(); it ++ ,i++ ){
	RedBaryCenConM::computeConM(*it,U,con_Ms[i]);
  }
}

