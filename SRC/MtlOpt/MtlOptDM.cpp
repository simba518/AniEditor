#include "MtlOptDM.h"
using namespace MTLOPT;

MtlOptDM::MtlOptDM(){
  T = 0;
  T_begin = 0;
  T_end = 0;
  alphaK0 = 0.0f;
  alphaM0 = 0.0f;
  mechanic_con_frames = 2;
}

void MtlOptDM::setTotalFrames(const int T){
  assert_ge(T,3);
  this->T = T;
  T_end = T-1;
}
void MtlOptDM::setDamping(const double ak,const double am){
  assert_ge(ak,0);
  assert_ge(am,0);
  alphaK0 = ak;
  alphaM0 = am;
}
void MtlOptDM::initVariables(const VectorXd &lambda0, const int rs){
	  
  assert_gt(rs,0);
  assert_ge(lambda0.size(),rs);
  assert_ge(T,3);

  Lambda0 = lambda0;
  Lambda = lambda0.head(rs);
  alphaM = alphaM0*VectorXd::Ones(rs);
  alphaK = alphaK0*VectorXd::Ones(rs);
  Damping = alphaM+alphaK.asDiagonal()*Lambda;

  Z.resize(rs, T);
  Z.setZero();
  S = MatrixXd::Identity(lambda0.size(), rs);
  K = Lambda.asDiagonal();
}

void MtlOptDM::fixHeadFrames(const int numFrames){
  assert_ge(numFrames,0);
  T_begin = numFrames;
}
void MtlOptDM::fixTailFrames(const int numFrames){
  assert_ge(numFrames,0);
  assert_lt(numFrames,T);
  T_end = T-numFrames-1;
}
void MtlOptDM::setMechanicConFrames(const int numFrames){
  assert_in(numFrames,0,subFrames()-1);
  mechanic_con_frames = numFrames;
}

void MtlOptDM::setConGroups(const int frame_id,const vector<int>&g,const VectorXd&pc){

  assert_in(frame_id,T_begin,T_end);
  assert_eq(g.size()*3, pc.size());

  removePartialCon(frame_id);

  conFrames.push_back(frame_id);
  conNodes.push_back(g);
  uc.push_back(pc);
  sortConGroups();
}
void MtlOptDM::setUc(const int frame_id, const VectorXd &pc){
	  
  assert_in(frame_id,T_begin,T_end);
  assert_eq(pc.size()%3, 0);
  int i = 0;
  for (i = 0; i < (int)conFrames.size(); ++i){
	if (frame_id == conFrames[i]){
	  assert_eq(this->uc[i].size(),pc.size());
	  this->uc[i] = pc;
	  break;
	}
  }
  ERROR_LOG_COND("failed to setUc for frame "<< frame_id, (i < (int)conFrames.size()));
}
void MtlOptDM::removeAllPosCon(){
  conNodes.clear();
  conFrames.clear();
  uc.clear();
}
void MtlOptDM::removePartialCon(const int frame_id){

  const vector<int> conFrames_t = conFrames;
  const vector<vector<int> > conNodes_t = conNodes;
  const vector<VectorXd> uc_t = uc;

  conFrames.clear();
  conNodes.clear();
  uc.clear();
  for (size_t i = 0; i < conFrames_t.size(); ++i){
	if (frame_id != conFrames_t[i]){
	  conFrames.push_back(conFrames_t[i]);
	  conNodes.push_back(conNodes_t[i]);
	  uc.push_back(uc_t[i]);
	}
  }
}
void MtlOptDM::sortConGroups(){

  set<std::pair<int,int> > con;
  for (int i = 0; i < conFrames.size(); ++i){
	con.insert(std::pair<int,int>(conFrames[i],i));
  }

  const vector<vector<int> > conNodes_t = conNodes;
  const vector<VectorXd> uc_t = uc;

  uc.clear();
  conNodes.clear();
  conFrames.clear();
  typedef std::pair<int,int> value_type; 
  BOOST_FOREACH(const value_type &p, con){
	const int i = p.second;
	const int f = p.first;
	uc.push_back(uc_t[i]);
	conNodes.push_back(conNodes_t[i]);
	conFrames.push_back(f);
  }
}

void MtlOptDM::setKeyfames(const vector<int> &keyf, const MatrixXd &Zzz){

  assert_eq(keyf.size(),Zzz.cols());
  assert_eq(Zzz.rows(),S.rows());
  set<std::pair<int,int> > con;
  for (int i = 0; i < keyf.size(); ++i){
	con.insert(std::pair<int,int>(keyf[i],i));
  }
  assert_eq(con.size(),keyf.size());

  Zk = Zzz;
  keyframes.clear();
  typedef std::pair<int,int> value_type; 
  BOOST_FOREACH(const value_type &p, con){
	const int i = p.second;
	const int f = p.first;
	Zk.col(keyframes.size()) = Zzz.col(i);
	keyframes.push_back(f);
	assert_gt(Z.cols(),f);
	assert_le(Z.rows(), Zzz.rows());
	// Z.col(f) = Zzz.col(i).head(Z.rows());
  }
}

int MtlOptDM::numConNodes()const{
  int n = 0;
  for (int i = 0; i < conNodes.size(); ++i)
	n += conNodes[i].size();
  return n;
}
int MtlOptDM::reducedDim()const{
  return Z.rows();
}
int MtlOptDM::subFrames()const{
  return (T_end-T_begin+1);
}

void MtlOptDM::updateZ(const double *subZ){

  assert_eq(Z.cols(),T);
  if (T_begin > 0){
	Z.leftCols(T_begin).setZero();
  }
  if (T_end < T-1){
	Z.rightCols(T-(T_end+1)).setZero();
  }
  assert(subZ != &Z(0,0)); // can not copy to it self.
  memcpy(&Z(0,T_begin), subZ, subFrames()*reducedDim()*sizeof(double));
}
int MtlOptDM::expandS(const int rs){

  if (rs == S.cols())
	return rs;

  
  // { // S = [S, O]
  // 	const int rw = S.rows();
  // 	assert_gt( rs, S.cols() );
  // 	assert_ge( rw, rs );
  // 	assert_gt( S.cols(), 0);
  // 	const MatrixXd ts = S;
  // 	S.resize(rw,rs);
  // 	S.setZero();
  // 	S.leftCols( ts.cols() ) = ts;  
  // }

  { // S = [S, I]
  	const int rw = S.rows();
  	assert_gt( rs, S.cols() );
  	assert_ge( rw, rs );
  	assert_gt( S.cols(), 0);
  	const MatrixXd ts = S;
  	S = MatrixXd::Identity(rw,rs);
  	S.leftCols( ts.cols() ) = ts;
  }
	  
  // Z = (Z^t, O)^t
  const MatrixXd tZ = Z;
  Z.resize( rs, Z.cols());
  Z.setZero();
  Z.topRows( tZ.rows() ) = tZ;

  // Lambda = [Lambda, Lambda0]
  assert_gt( Lambda.size(), 0 );
  assert_ge( Lambda0.size(), rs);
  const VectorXd tL = Lambda;
  Lambda.resize(rs);
  Lambda.head( tL.size() ) = tL;
  for (int i = tL.size(); i < rs; ++i){
	Lambda[i] = Lambda0[i];
  }
  K = Lambda.asDiagonal();

  // alphaM = alphaM0*VectorXd::Ones(rs);
  // alphaK = alphaK0*VectorXd::Ones(rs);
  // Damping = alphaM+alphaK.asDiagonal()*Lambda;

  return rs;
}
void MtlOptDM::decomposeK(){
	  
  assert_eq(K.rows(), K.cols());
  assert_eq(K.rows(), reducedDim());
  SelfAdjointEigenSolver<MatrixXd> eigenK(K);
  const MatrixXd &U = eigenK.eigenvectors();
  Z = U.transpose()*Z;
  S = S*U;
  Lambda = eigenK.eigenvalues();
  K = Lambda.asDiagonal();
}

void MtlOptDM::print()const{

  INFO_LOG("rs: " << S.cols());
  INFO_LOG("rw: " << S.rows());
  INFO_LOG("old_alpha_k: " << alphaK0);
  INFO_LOG("old_alpha_m: " << alphaM0);
  INFO_LOG("new_damping: " << Damping.transpose());
  INFO_LOG("new_lambda: " << Lambda.transpose());
  // INFO_LOG("K: " << K);
  INFO_LOG("T_begin: " << T_begin);
  INFO_LOG("T_end: " << T_end);
  INFO_LOG("num_of_con_frames: " << conFrames.size());
  INFO_LOG("num_of_con_nodes: " << numConNodes());
}
