#include <MatrixIO.h>
#include <MatrixTools.h>
#include <volumetricMeshLoader.h>
#include <RedCenConM.h>
#include <FullBaryCenConM.h>
#include <JsonFilePaser.h>
#include <Log.h>
#include "RedRSInterpolator.h"
using namespace UTILITY;
using namespace LSW_ANI_EDITOR;

RedRSInterpolator::RedRSInterpolator(){

  warper = pRSWarperExt(new RSWarperExt());
  use_warp = true;
}

bool RedRSInterpolator::init (const string init_filename){
  
  TRACE_FUN();
  bool succ = anieditor.initialize(init_filename);
  if(!succ){
	ERROR_LOG("failed to initialize RedRSInterpolator::anieditor.");
  } else{

	JsonFilePaser json_f;
	succ = json_f.open(init_filename);
	if (succ){

	  // read linear basis 
	  succ &= json_f.readMatFile("eigen_vectors",W);

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
	}
	if(succ){
	  anieditor.setTotalFrame(u_ref.size());
	}

	//@todo
	succ &= json_f.readMatFile("nonlinear_basis",B);

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

  if (!succ){
	ERROR_LOG("failed to initialize RedRSInterpolator.");
  }
  return succ;
}

void RedRSInterpolator::setConGroups(const int frame_id,const vector<set<int> >&group, const VectorXd &uc){

  addConGroups(frame_id,group,uc);
  anieditor.setConstrainedFrames(con_frame_id);
}

void RedRSInterpolator::setUc(const int frame_id, const VectorXd &uc){

  assert_in (frame_id,0,(int)u_ref.size()-1);
  assert_eq (con_nodes.size(),con_frame_id.size());
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
  ERROR_LOG_COND("failed to setUc!", (i < (int)con_frame_id.size())); 
}

void RedRSInterpolator::removeAllPosCon(){
  
  con_frame_id.clear();
  con_nodes.clear();
  uc.clear();
  /// @todo remove constraints in anieditor.
}

bool RedRSInterpolator::interpolate (){

  if (modalDisplayer.isShowModalModes()){
	return modalDisplayer.showModalModes(anieditor.getEigenvalues(),delta_z);
  }

  vector<MatrixXd > Cs(con_frame_id.size());
  const int r = reducedDim();
  const int iterNum = 6; ///@todo
  if (con_frame_id.size() > 0){
	for (int it = 0; it < iterNum; ++it){

	  // get partial constraints
	  VectorXd new_uc;
	  for (size_t f = 0; f < con_frame_id.size(); ++f){

		Cs[f].resize(con_nodes[f].size()*3,r);
		const VectorXd &z = delta_z[ con_frame_id[f] ];
		for (size_t i = 0; i < con_nodes[f].size(); ++i){

		  const int con_node = con_nodes[f][i];
		  const Vector3d u3 = nodeWarper->warp(z,con_frame_id[f],con_node);
		  Eigen::Matrix<double,3,-1> C;
		  nodeWarper->jacobian(z,f,con_node,C);
		  const Vector3d c3 = uc[f].segment(i*3,3)+C*z-u3;
		  new_uc.conservativeResize(new_uc.size()+3);
		  new_uc.segment(new_uc.size()-3,3) = c3;
		  Cs[f].block(3*i,0,3,r) = C;
		}
	  }

	  // solve
	  VectorXd zc;
	  anieditor.computeConstrainedFrames(Cs,new_uc,zc);

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

const vector<VectorXd>& RedRSInterpolator::getDeltaInterpSeq()const{
  return delta_z;
}

bool RedRSInterpolator::validConstraints(const int frame_id)const{

  const bool valid = !( anieditor.isBoundaryFrame(frame_id) );
  WARN_LOG_COND("the frame id of the constraints is invalid:"<<frame_id,valid);
  return valid;
}

const VectorXd& RedRSInterpolator::getInterpU(const int frame_id){

  assert_in (frame_id, 0, getT()-1);
  if (use_warp && warper != NULL){
	nodeWarper->warp(delta_z[frame_id],frame_id,full_u);
	// VectorXd tempt_u;
	// const VectorXd p = W*delta_z[frame_id];
	// warper->warp(p,frame_id,tempt_u);
	// const double vol = tet_mesh->getTotalVolume();
	// cout << "RS diff REDUCED-RS: " << (full_u-tempt_u).norm()/vol << endl;
  }else{
	// const VectorXd p = W*delta_z[frame_id];
	// warper->warp(p,frame_id,full_u);
	full_u = u_ref[frame_id] + W*delta_z[frame_id];
  }
  return full_u;
}

const VectorXd& RedRSInterpolator::getInputU(const int frame_id){

  assert_in (frame_id, 0, getT()-1);
  return u_ref[frame_id];
}

const VectorXd& RedRSInterpolator::getUnwarpedU(const int frame_id){

  unwarped_full_u = u_ref[frame_id] + W * delta_z[frame_id];
  return unwarped_full_u;
}

bool RedRSInterpolator::initWarper(const string init_filename){

  JsonFilePaser inf;
  bool succ = inf.open(init_filename);
  if (succ){
	
  	// load and set volumetric mesh
  	string vol_filename;
  	succ &= inf.readFilePath("vol_filename",vol_filename);

	tet_mesh = pVolumetricMesh(VolumetricMeshLoader::load(vol_filename.c_str()));
  	succ &= (tet_mesh!=NULL);

	if (tet_mesh != NULL){
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
  string cubPstr,cubWstr;
  vector<int> cubP;
  vector<double> cubW;
  if (inf.readFilePath("cubaturePoints",cubPstr)&&
	  inf.readFilePath("cubatureWeights",cubWstr) &&
	  UTILITY::loadVec(cubPstr,cubP,UTILITY::TEXT)&&
	  UTILITY::loadVec(cubWstr,cubW)){
	nodeWarper = pRedRSWarperExt(new RedRSWarperExt(warper,B,W,cubP,cubW));
  }else{
	nodeWarper = pRedRSWarperExt(new RedRSWarperExt(warper,B,W));
  }
  return succ;
}

bool RedRSInterpolator::editable(const int frame_id)const{
  
  return !(anieditor.isBoundaryFrame(frame_id));
}

const VectorXd &RedRSInterpolator::getWarpU(const int frame_id,
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

void RedRSInterpolator::setAllConGroups(const set<pConNodesOfFrame> &newCons){

  removeAllPosCon();
  BOOST_FOREACH(pConNodesOfFrame con, newCons){
  	if( !con->isEmpty() ){
	  addConGroups(con->getFrameId(), con->getConNodesSet(), con->getBarycenterUc());
	}
  }
  anieditor.setConstrainedFrames(con_frame_id);
}

void RedRSInterpolator::addConGroups(const int frame_id,const vector<set<int> >&group,const VectorXd&uc){
  
  if(!validConstraints(frame_id)){
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

void RedRSInterpolator::print(const string data_name)const{
  
  if(string("eigenvalues") == data_name){
	cout<< "eigenvalues: \n" << anieditor.getEigenvalues() << endl;
  }else{
	BaseInterpolator::print(data_name);
  }
}

void RedRSInterpolator::removeConOfFrame(const int frame_id){
  
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
