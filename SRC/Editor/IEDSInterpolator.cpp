#include <MatrixIO.h>
#include <MatrixTools.h>
#include <ConMatrixTools.h>
#include <JsonFilePaser.h>
#include <RSWarperExt.h>
#include <ModalWarpingExt.h>
#include <RotationMatrix.h>
using namespace UTILITY;
using namespace LSW_WARPING;

#include "IEDSInterpolator.h"
using namespace LSW_ANI_EDITOR;

IEDSInterpolator::IEDSInterpolator(){

  warper = pRSWarperExt(new RSWarperExt());
  use_warp = true;
}

bool IEDSInterpolator::init (const string init_filename){
  
  TRACE_FUN();
  bool succ = anieditor.initialize(init_filename);
  if(!succ){
	ERROR_LOG("failed to initialize IEDSInterpolator::anieditor.");
  } else{

	JsonFilePaser json_f;
	succ = json_f.open(init_filename);
	if (succ){
	  // read linear basis 
	  succ = json_f.readMatFile("eigen_vectors",W);
	}

	// read input sequence
	int zero_input_animation = 0;
	if (json_f.read("zero_input_animation",zero_input_animation) && zero_input_animation > 0){
	  u_ref.resize(zero_input_animation);
	  for (size_t i = 0; i < u_ref.size(); ++i){
		u_ref[i].resize(W.rows());
		u_ref[i].setZero();
	  }
	}else{
	  succ = false;
	  string fullinput;
	  if(json_f.readFilePath("full_input_animation",fullinput)){
		MatrixXd Uin;
		succ = EIGEN3EXT::load(fullinput,Uin);
		EIGEN3EXT::convert(Uin,u_ref);
	  }
	}

	// read in refrence sequence
	if(succ){
	  anieditor.setTotalFrame(u_ref.size());
	}
	  
	// read in nonlinear and linear basis 
	if(succ){
	  succ &= json_f.readMatFile("eigen_vectors",W);
	}

	// initialize delta_z to be zeros.
	const int T = this->getT();
	const int decoupled_r = anieditor.reducedDim();
	delta_z.resize(T);
	for (int i = 0; i < T; ++i){
	  delta_z[i].resize(decoupled_r);
	  delta_z[i].setZero();
	}
  }
  if (succ){
	succ = anieditor.precompute();
	succ &= initWarper(init_filename);
  }
  ERROR_LOG_COND("failed to initialize IEDSInterpolator.",succ);
  return succ;
}

void IEDSInterpolator::setConGroups(const int frame_id,const vector<set<int> >&group, const VectorXd &uc){

  if(!validConstraints(frame_id)){
	return ;
  }
  assert_in (frame_id,0,getT()-1);

  if (group.size() <= 0){
	removeConOfFrame(frame_id);
	return ;
  }

  SparseMatrix<double> C;
  computeBaryCenterConM(group,W.rows()/3, C);

  int i = 0;
  for (i = 0; i < (int)con_frame_id.size(); ++i){

	if (frame_id == con_frame_id[i]){
	  // frame_id is in con_frame_id
	  this->C[i] = C;
	  this->uc[i] = uc;
	  break;
	}
  }

  // frame_id is not in con_frame_id
  if ( i >= (int)con_frame_id.size()){

	con_frame_id.push_back(frame_id);
	this->C.push_back(C);
	this->uc.push_back(uc);
  }
}

void IEDSInterpolator::setUc(const int frame_id, const VectorXd &uc){

  assert_in (frame_id,0,(int)u_ref.size()-1);
  assert_eq (C.size(),con_frame_id.size());
  if(!validConstraints(frame_id)){
	return ;
  }

  int i = 0;
  for (i = 0; i < (int)con_frame_id.size(); ++i) {

	if (frame_id == con_frame_id[i]){
	  this->uc[i] = uc;
	  break;
	}
  }
  ERROR_LOG_COND("failed to setUc!", (i<(int)con_frame_id.size())); 
}

void IEDSInterpolator::removeAllPosCon(){
  
  con_frame_id.clear();
  C.clear();
  uc.clear();
}

bool IEDSInterpolator::interpolate (){

  // set partial constraints
  applyPosCon();

  // solve
  bool succ = anieditor.applyEdits();

  // get results
  if (succ){
  	const int T = getT();
  	const int decoupled_r = anieditor.reducedDim();
  	const VectorXd &Z = anieditor.getEditResult();
  	assert_eq(Z.size(),T*decoupled_r);
  	delta_z.resize(T);
  	for (int i = 0; i < T; ++i){
  	  delta_z[i] = Z.segment(i*decoupled_r, decoupled_r);
  	}
  }

  return succ;
}

const vector<VectorXd>& IEDSInterpolator::getDeltaInterpSeq()const{
  
  return delta_z;
}

bool IEDSInterpolator::validConstraints(const int frame_id)const{

  const bool valid = !( anieditor.isRemoved(frame_id) );
  WARN_LOG_COND("the frame id of the constraints is invalid:"<<frame_id, valid);
  return valid;
}

const VectorXd& IEDSInterpolator::getInterpU(const int frame_id){

  assert_in (frame_id, 0, getT()-1);
  if (use_warp && warper != NULL){
	const VectorXd p = W*delta_z[frame_id];
	warper->warp(W*delta_z[frame_id], frame_id,full_u);
  }else{
	full_u = u_ref[frame_id] + W*delta_z[frame_id];
  }
  return full_u;
}

const VectorXd& IEDSInterpolator::getInputU(const int frame_id){

  assert_in (frame_id, 0, getT()-1);
  return u_ref[frame_id];
}

const VectorXd& IEDSInterpolator::getUnwarpedU(const int frame_id){

  unwarped_full_u = u_ref[frame_id] + W * delta_z[frame_id];
  return unwarped_full_u;
}

bool IEDSInterpolator::initWarper(const string init_filename){
  
  JsonFilePaser inf;
  bool succ = inf.open(init_filename);
  if (succ){
	
  	// load and set volumetric mesh
  	string vol_filename;
  	succ &= inf.readFilePath("vol_filename",vol_filename);

  	pTetMesh tet_mesh = pTetMesh(new TetMesh());
	succ = tet_mesh->load(vol_filename);
	if (succ){
	  warper->setTetMesh(tet_mesh);
	}

  	// load and set fixed nodes
  	string fixed_nodes_RS2Euler;
  	if( inf.readFilePath("fixed_nodes_RS2Euler",fixed_nodes_RS2Euler) ){
	  vector<int> fixed_nodes;
	  if(UTILITY::loadVec(fixed_nodes_RS2Euler,fixed_nodes,UTILITY::TEXT))
		warper->setFixedNodes(fixed_nodes);
	}

  	// set reference displacements
	warper->setRefU(u_ref);

  	// precompute
	succ &= warper->precompute();

	// use warp
	if ( !inf.read("use_warp",use_warp) ){
	  use_warp = false;
	}
  }
  if (!succ) {
	use_warp = false;
  }
  return succ;
}

bool IEDSInterpolator::editable(const int frame_id)const{
  
  return anieditor.editable(frame_id);
}

const VectorXd &IEDSInterpolator::getWarpU(const int frame_id,
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

void IEDSInterpolator::applyPosCon(){

  // method 1
  vector<MatrixXd > Cw; 
  vector<VectorXd> delta_uc; 
  for (size_t i = 0; i < con_frame_id.size(); ++i){

  	Cw.push_back(C[i]*W);
  	delta_uc.push_back(uc[i] - C[i]*u_ref[con_frame_id[i]]);
  }
  anieditor.setPosCon(con_frame_id,Cw, delta_uc);
}

void IEDSInterpolator::removeConOfFrame(const int frame_id){
  
  assert_in(frame_id,0,getT()-1);
  assert_eq(con_frame_id.size(), C.size());
  assert_eq(con_frame_id.size(), uc.size());

  const vector<int> con_frame_id_t = con_frame_id;
  const vector<SparseMatrix<double> > C_t = C;
  const vector<VectorXd> uc_t = uc; 
  con_frame_id.clear();
  C.clear();
  uc.clear();

  for (size_t i = 0; i < con_frame_id_t.size(); ++i){
	if (frame_id != con_frame_id_t[i]){
	  con_frame_id.push_back(con_frame_id_t[i]);
	  C.push_back(C_t[i]);
	  uc.push_back(uc_t[i]);
	}
  }
}
