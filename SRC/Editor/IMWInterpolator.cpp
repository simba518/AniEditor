#include <ConMatrixTools.h>
#include <JsonFilePaser.h>
#include <RotationMatrix.h>
#include "IMWInterpolator.h"
#include <MatrixIO.h>
#include <MatrixTools.h>
using namespace UTILITY;
using namespace LSW_WARPING;
using namespace LSW_ANI_EDITOR;

IMWInterpolator::IMWInterpolator(){

  warper = pModalWarpingExt(new ModalWarpingExt());
  use_warp = true;
}

bool IMWInterpolator::init (const string init_filename){
  
  TRACE_FUN();

  bool succ = anieditor.initialize(init_filename);
  if(!succ){
	ERROR_LOG("failed to initialize IMWInterpolator::anieditor.");
  } else{

	JsonFilePaser json_f;
	succ = json_f.open(init_filename);
	if (succ){
	  // read linear basis 
	  succ &= json_f.readMatFile("eigenvectors",W);
	  loadUref(json_f,u_ref);
	}
	if (succ){
	  anieditor.setTotalFrame(u_ref.size());
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
	anieditor.precompute();
	succ = initWarper(init_filename);
  }

  modalDisplayer.initialize(init_filename);
  ERROR_LOG_COND("failed to initialize IMWInterpolator.",succ);
  return succ;
}

void IMWInterpolator::setConGroups(const int frame_id,const vector<set<int> >&group, const Matrix<double,3,-1> &uc){

  vector<set<int> > group_rlst;
  VectorXd uc_rlst;
  splitAllConstraints(group,uc,group_rlst,uc_rlst);
  addConGroups(frame_id,group_rlst,uc_rlst);
  anieditor.setConstrainedFrames(con_frame_id);
}
 
void IMWInterpolator::setUc(const int frame_id, const Matrix<double, 3,-1> &uc){

  assert_in (frame_id,0,(int)u_ref.size()-1);
  assert_eq (con_nodes.size(),con_frame_id.size());
  if(!validConstraints(frame_id)){
	return ;
  }

  int i = 0;
  for (i = 0; i < (int)con_frame_id.size(); ++i) {
	if (frame_id == con_frame_id[i]){
	  assert_eq(this->uc[i].size(),uc.size());
	  this->uc[i] = Map<VectorXd>(const_cast<double*>(&uc(0,0)),uc.size());
	  break;
	}
  }
  ERROR_LOG_COND("failed to setUc!",(i < (int)con_frame_id.size())); 
}

void IMWInterpolator::removeAllPosCon(){
  
  con_frame_id.clear();
  con_nodes.clear();
  uc.clear();
  /// @todo remove constraints in anieditor.
}

bool IMWInterpolator::interpolate (){

  if (modalDisplayer.isShowModalModes()){
	return modalDisplayer.showModalModes(anieditor.getEigenvalues(),delta_z);
  }

  vector<MatrixXd > Cs(con_frame_id.size());
  Matrix3d Jw,Jp;
  Vector3d u3;
  const int r = reducedDim();
  const int iterNum = 5;
  if (con_frame_id.size() > 0){

	for (int it = 0; it < iterNum; ++it){

	  // get partial constraints
	  VectorXd new_uc;
	  for (size_t f = 0; f < con_frame_id.size(); ++f){

		Cs[f].resize(con_nodes[f].size()*3,r);
		const VectorXd &z = delta_z[ con_frame_id[f] ];
		const VectorV3 &inputW = warper->getInputRotVec(con_frame_id[f]);
		const VectorXd &inputQ = warper->getInputLocalDisp(con_frame_id[f]);
		for (size_t i = 0; i < con_nodes[f].size(); ++i){

		  const int con_node = con_nodes[f][i];
		  Vector3d w3 = inputW[con_node] + rotW.block(con_node*3,0,3,r)*z;
		  if (w3.norm() <= 1e-8){
			w3[0] = 1e-8;
		  }
		  const Vector3d p3 = inputQ.segment(con_node*3,3) + W.block(con_node*3,0,3,r)*z;
		  NodeWarper::warp(w3,p3,u3);
		  warperExtAD.jacobian(w3,p3,con_node,Jw,Jp);
		  const Matrix<double,3,-1> C = Jw*rotW.block(con_node*3,0,3,r) + Jp*W.block(con_node*3,0,3,r);
		  const Vector3d c3 = uc[f].segment(i*3,3)+C*z-u3;
		  new_uc.conservativeResize(new_uc.size()+3);
		  new_uc.segment(new_uc.size()-3,3) = c3;
		  Cs[f].block(3*i,0,3,r) = C;
		}
	  }

	  // solve
	  VectorXd zc;
	  Timer timer;
	  timer.start();
	  anieditor.computeConstrainedFrames(Cs,new_uc,zc);
	  timer.stop("KKT: ");

	  assert_eq(zc.size(), con_frame_id.size()*r);
	  for (size_t i = 0; i < con_frame_id.size(); ++i){
		delta_z[con_frame_id[i]] = zc.segment(r*i,r);
	  }

	  /// @todo test
	  static VectorXd lastZ0;
	  assert_gt (con_frame_id.size(),0);
	  if ( lastZ0.size() == delta_z[con_frame_id[0]].size() ){
		cout<< "z diff: " << (lastZ0 - delta_z[con_frame_id[0]]).norm() << endl;
	  }
	  lastZ0 = delta_z[con_frame_id[0]];
	}
  }

  // get the results
  anieditor.updateZ();
  const int T = getT();
  const VectorXd &z = anieditor.getEditResult();
  assert_eq(z.size(),(T-anieditor.numBoundaryFrames())*r);
  delta_z.resize(T);
  int count = 0;
  for (int i = 0; i < T; ++i){
  	if ( !anieditor.isBoundaryFrame(i) ){
  	  delta_z[i] = z.segment(count*r,r);
  	  count ++;
  	}else{
  	  delta_z[i].setZero();
  	}
  }
  return true;
}

const vector<VectorXd>& IMWInterpolator::getDeltaInterpSeq()const{
  return delta_z;
}

bool IMWInterpolator::validConstraints(const int frame_id)const{

  const bool valid = !( anieditor.isBoundaryFrame(frame_id) );
  WARN_LOG_COND("the frame id of the constraints is invalid:"<<frame_id,valid);
  return valid;
}

const VectorXd& IMWInterpolator::getInterpU(const int frame_id){

  assert_in (frame_id, 0, getT()-1);
  if (use_warp && warper != NULL){
	Timer timer;
	timer.start();
	const VectorXd p = W*delta_z[frame_id];
	timer.stop("p=Wz: ");
	timer.start();
	warper->warp(p,frame_id,full_u);
	timer.stop("MW: ");
  }else{
	full_u = u_ref[frame_id] + W*delta_z[frame_id];
  }
  return full_u;
}

const VectorXd& IMWInterpolator::getInputU(const int frame_id){

  assert_in (frame_id, 0, getT()-1);
  return u_ref[frame_id];
}

const VectorXd& IMWInterpolator::getUnwarpedU(const int frame_id){

  unwarped_full_u = u_ref[frame_id] + W * delta_z[frame_id];
  return unwarped_full_u;
}

bool IMWInterpolator::initWarper(const string init_filename){

  TRACE_FUN();
  JsonFilePaser inf;
  bool succ = inf.open(init_filename);
  pTetMesh tet_mesh = pTetMesh(new TetMesh());
  if (succ){
	
  	// load and set volumetric mesh
  	string vol_filename;
 	succ &= inf.readFilePath("vol_file",vol_filename);

	succ = tet_mesh->load(vol_filename);
	if (succ){
	  warper->setTetMesh(tet_mesh);
	}

  	// set reference displacements
	warper->setRefU(u_ref);

  	// precompute
	Timer timer;
	timer.start();
	succ &= warper->precompute();
	timer.stop("initialize warper: ");

	// use warp
	if ( !inf.read("use_warp",use_warp) ){
	  use_warp = false;
	}
  }
  if (!succ){
	use_warp = false;
  }

  Timer timer;
  timer.start();
  warperExtAD.precompute(tet_mesh);
  timer.stop("initialize warperext: ");

  if (!inf.readMatFile("rot_modal_matrix",rotW)){
	INFO_LOG("compute rotational modal warping matrix....");
	timer.start();
	WarperAD::computeRotMat(tet_mesh,W,rotW);
	timer.stop("initialize rot modal matrix: ");
  }
  assert_eq(rotW.rows(),W.rows());
  assert_eq(rotW.cols(),W.cols());
  return succ;
}

bool IMWInterpolator::editable(const int frame_id)const{
  
  return !(anieditor.isBoundaryFrame(frame_id));
}

const VectorXd &IMWInterpolator::getWarpU(const int frame_id,
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

void IMWInterpolator::addConGroups(const int frame_id,const vector<set<int> >&group,const VectorXd&uc){

  if(!validConstraints(frame_id)){
	return ;
  }
  assert_in (frame_id,0,getT()-1);
  assert_eq (uc.size()%3,0);
  assert_eq (uc.size(),group.size()*3);

  if(group.size() <= 0){
	removeConOfFrame(frame_id);
	return ;
  }

  vector<int> con_nodes;
  for (size_t i = 0; i < group.size(); ++i){
	assert_gt(group[i].size() , 0);
	con_nodes.push_back(*(group[i].begin()));
  }

  // we need to sort the constraints by frames numder, required by anieditor.
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

void IMWInterpolator::setConGroups(const int frame_id,const vector<set<int> >&group,const VectorXd&uc){
  
  if(!validConstraints(frame_id)){
	return ;
  }
  assert_in (frame_id,0,getT()-1);
  assert_eq (uc.size()%3,0);
  assert_eq (uc.size(),group.size()*3);

  if(group.size() <= 0){
	removeConOfFrame(frame_id);
	return ;
  }

  vector<int> con_nodes;
  for (size_t i = 0; i < group.size(); ++i){
	assert_gt(group[i].size() , 0);
	con_nodes.push_back(*(group[i].begin()));
  }

  // we need to sort the constraints by frames numder, required by anieditor.
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

bool IMWInterpolator::save(const string data_name, const string file_name)const{
  
  if (string("rotational_modal_matrix") == data_name){
	return EIGEN3EXT::write(file_name,rotW);
  }
  return BaseInterpolator::save(data_name,file_name);
}

void IMWInterpolator::print(const string data_name)const{
  
  if(string("eigenvalues") == data_name){
	cout<< "eigenvalues: \n" << anieditor.getEigenvalues() << endl;
  }else{
	BaseInterpolator::print(data_name);
  }
}

void IMWInterpolator::removeConOfFrame(const int frame_id){
  
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

bool IMWInterpolator::loadUref(JsonFilePaser &json_f,vector<VectorXd> &u_ref)const{

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
