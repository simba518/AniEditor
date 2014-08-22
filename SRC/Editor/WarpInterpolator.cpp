#include <MatrixIO.h>
#include <MatrixTools.h>
#include <ConMatrixTools.h>
#include <JsonFilePaser.h>
#include <Log.h>
#include "WarpInterpolator.h"
using namespace UTILITY;
using namespace LSW_ANI_EDITOR;

bool WarpInterpolator::init (const string init_filename){
  
  TRACE_FUN();
  JsonFilePaser json_f;
  VectorXd lambda;
  bool succ = json_f.open(init_filename);
  if (succ){

	succ &= json_f.readVecFile("eigenvalues",lambda);
	succ &= json_f.readMatFile("eigenvectors",W);
	succ &= json_f.readMatFile("nonlinear_basis",B);
	assert_eq(lambda.size(),W.cols());

	int use_sub_MA_modes = 0;
	if (json_f.read("rw",use_sub_MA_modes) && use_sub_MA_modes > 0){
	  assert_in(use_sub_MA_modes,1,lambda.size());
	  const MatrixXd tW = W;
	  W = tW.leftCols(use_sub_MA_modes);
	  Lambda = lambda.head(use_sub_MA_modes);
	}else{
	  Lambda = lambda;
	}
	use_warp = json_f.read("use_warp",use_warp) ? use_warp:false;
	succ &= loadUref(json_f,u_ref);
  }
  if (succ){
	delta_z.resize(getT());
	for (int i = 0; i < getT(); ++i){
	  delta_z[i].resize(linearSubDim());
	  delta_z[i].setZero();
	}
  }
  if (succ){
	succ = initWarper(json_f);
  }
  modalDisplayer.initialize(init_filename);
  if (!succ) ERROR_LOG("failed to initialize WarpInterpolator.");
  return succ;
}

void WarpInterpolator::setUc(const int frame_id, const Matrix<double,3,-1>&uc){

  assert_in (frame_id,0,(int)u_ref.size()-1);
  assert_eq (con_nodes.size(),con_frame_id.size());
  if(!validConstraints(frame_id)){
	ERROR_LOG("frame "<<frame_id<<" can not constrained.");
	return ;
  }
  int i = 0;
  for (i = 0; i < (int)con_frame_id.size(); ++i){
	if (frame_id == con_frame_id[i]){
	  assert_eq(this->uc[i].size(),uc.size());
	  assert_ge(uc.size(),3);
	  this->uc[i] = Map<VectorXd>(const_cast<double*>(&uc(0,0)),uc.size());
	  break;
	}
  }
  ERROR_LOG_COND("failed to setUc for frame "<< frame_id, (i < (int)con_frame_id.size()));
}

const VectorXd& WarpInterpolator::getInterpU(const int frame_id){

  assert_in (frame_id, 0, getT()-1);
  if (use_warp && nodeWarper != NULL){
	nodeWarper->warp(delta_z[frame_id],frame_id,full_u);
  }else{
	const VectorXd p = W*delta_z[frame_id];
	// full_u = u_ref[frame_id] + p;
	warper->warp(p,frame_id,full_u); /// @todo
  }
  return full_u;
}

void WarpInterpolator::sortCubPoints(vector<int> &cubPoints, vector<double> &cubWeights){

  assert_eq(cubPoints.size(), cubWeights.size());
  set<pair<int,double> > s;
  for (int i = 0; i < cubPoints.size(); ++i){
    s.insert(pair<int,double>(cubPoints[i], cubWeights[i]));
  }

  cubPoints.clear();
  cubWeights.clear();
  set<pair<int,double> >::iterator it = s.begin();
  for (; it != s.end(); ++it){
	DEBUG_LOG(it->first <<", " << it->second);
    cubPoints.push_back(it->first);
	cubWeights.push_back(it->second);
  }
}

bool WarpInterpolator::initWarper(JsonFilePaser &inf){

  string vol_filename;
  bool succ = inf.readFilePath("vol_file",vol_filename);
  pTetMesh tet_mesh = pTetMesh(new TetMesh());
  if (succ) succ = tet_mesh->load(vol_filename);
  if (succ){ warper->setTetMesh(tet_mesh); }

  vector<int> fixed_nodes;
  if(inf.readVecFile("fixed_nodes",fixed_nodes,UTILITY::TEXT))
	warper->setFixedNodes(fixed_nodes);
  warper->setRefU(u_ref);
  if(succ) succ = warper->precompute();

  vector<int> cubP;
  vector<double> cubW;
  if (inf.readVecFile("cubature_points",cubP,UTILITY::TEXT)&&
	  inf.readVecFile("cubature_weights",cubW)){
	sortCubPoints(cubP, cubW);
  	nodeWarper = pRedRSWarperExt(new RedRSWarperExt(warper,B,W,cubP,cubW));
  }else{
  	nodeWarper = pRedRSWarperExt(new RedRSWarperExt(warper,B,W));
  }
  return succ;
}

const VectorXd &WarpInterpolator::getWarpU(const int frame_id,
											 const vector<set<int> > &c_nodes,
											 const VectorXd &bcen_uc){
  assert_in (frame_id, 0, getT()-1);
  assert (warper != NULL);
  warper->setConNodes(c_nodes,bcen_uc);
  warper->warp(W*delta_z[frame_id], frame_id,full_u);
  vector<set<int> > c_nodes_empty;
  VectorXd bcen_uc_emtpy;
  warper->setConNodes(c_nodes_empty,bcen_uc_emtpy);
  return full_u;
}

void WarpInterpolator::addConGroups(const int frame_id,const vector<set<int> >&group,const VectorXd &uc){

  if(!validConstraints(frame_id)){
	ERROR_LOG("frame "<<frame_id<<" can not be constrained.");
	return ;
  }
  assert_in (frame_id,0,getT()-1);
  assert_eq (uc.size()%3,0);
  assert_eq (uc.size(),group.size()*3);

  vector<int> con_nodes;
  for (size_t i = 0; i < group.size(); ++i){
	assert_gt(group[i].size() , 0);
	con_nodes.push_back(*(group[i].begin()));
  }

  if (con_nodes.size() <=0 ){
	removeConOfFrame(frame_id);
	return ;
  }

  // sort the constraints by frames numder.
  size_t i = 0;
  for (; i < con_frame_id.size(); ++i){
	if (frame_id == con_frame_id[i]){
	  // frame_id is in con_frame_id
	  this->con_nodes[i] = con_nodes;
	  this->uc[i] = uc;
	  break;
	}else if(frame_id < con_frame_id[i]){
	  this->con_nodes.insert(this->con_nodes.begin()+i,con_nodes);
	  this->uc.insert(this->uc.begin()+i,uc);
	  con_frame_id.insert(con_frame_id.begin()+i,frame_id);
	  break;
	}
  }
  if ( i >= con_frame_id.size()){
	this->con_nodes.push_back(con_nodes);
	this->uc.push_back(uc);
	con_frame_id.push_back(frame_id);
  }
}

void WarpInterpolator::removeConOfFrame(const int frame_id){
  
  assert_in(frame_id,0,getT()-1);
  assert_eq(con_frame_id.size(), con_nodes.size());
  assert_eq(con_nodes.size(),uc.size());
  
  const vector<int> con_frame_id_t = con_frame_id;   
  const vector<vector<int> > con_nodes_t = con_nodes;
  const vector<VectorXd> uc_t = uc; 
  con_frame_id.clear();
  con_nodes.clear();
  uc.clear();

  for (size_t i = 0; i < con_frame_id_t.size(); ++i){
	if (frame_id != con_frame_id_t[i]){
	  con_frame_id.push_back(con_frame_id_t[i]);
	  con_nodes.push_back(con_nodes_t[i]);
	  uc.push_back(uc_t[i]);
	}
  }
}

bool WarpInterpolator::loadUref(JsonFilePaser &json_f,vector<VectorXd> &u_ref)const{

  int T = 0;
  if( !json_f.read("T",T) )
	return false;
  assert_ge(T,3);

  MatrixXd W;
  if ( !json_f.readMatFile("eigenvectors",W) )
	return false;

  const int nx3 = W.rows();
  assert_gt(nx3,0);
  assert_eq(nx3%3,0);

  MatrixXd Uin;
  if ( json_f.readMatFile("input_animation",Uin) ){
	assert_eq(Uin.rows(), nx3);
	assert_ge(Uin.cols(), T);
	u_ref.resize(T);
	for (int i = 0; i < T; ++i) 
	  u_ref[i] = Uin.col(i);
  }else{
	u_ref.resize(T);
	for (int i = 0; i < T; ++i) {
	  u_ref[i].resize(nx3);
	  u_ref[i].setZero();
	}
  }

  return true;
}

VectorXd WarpInterpolator::adjustToSubdim(const VectorXd &z)const{
  assert_ge(z.size(),reducedDim());
  return z.head(reducedDim());
}

vector<VectorXd> WarpInterpolator::adjustToSubdim(const vector<VectorXd> &vz)const{
  vector<VectorXd> z = vz;
  for (int i = 0; i < z.size(); ++i)
	z[i] = adjustToSubdim(vz[i]);
  return z;
}
