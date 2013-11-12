#include <fstream>
#include <Timer.h>
#include <boost/foreach.hpp>
#include <MatrixIO.h>
#include <JsonFilePaser.h>
#include <InterpolatorFactory.h>
#include <MatrixTools.h>
#include "AniEditDM.h"
using namespace std;
using namespace UTILITY;
using namespace ANI_EDIT_UI;

AniEditDM::AniEditDM(pTetMeshEmbeding vol_obj, pAniDataModel animation):
  _animation(animation), vol_obj(vol_obj){

  record_drag = false;
}

AniEditDM::AniEditDM(pTetMeshEmbeding vol_obj, pAniDataModel animation, pBaseInterpolator interpolator):
  _animation(animation),interpolator(interpolator), vol_obj(vol_obj){

  record_drag = false;
}

bool AniEditDM::initialize(const string filename,const bool create_interp){

  TRACE_FUN();
  
  bool succ = false;
  if (interpolator == NULL || create_interp){
	interpolator = LSW_ANI_EDITOR::InterpolatorFactory::create(filename);
	succ = (interpolator != NULL);
  }else{
	succ = interpolator->init (filename);
  }

  if (succ){

	UTILITY::JsonFilePaser inf;
	if ( inf.open(filename) ){

	  string con_path, partial_con, drag_trajectory;
	  if ( inf.readFilePath("con_path", con_path) )
		succ = loadConPath(con_path);
	  if ( inf.readFilePath("partial_constraints",partial_con) )
		succ &= loadParitalCon(partial_con);
	  if ( inf.readFilePath("drag_trajectory",drag_trajectory) )
		succ &= loadDragRecord(drag_trajectory);
	  
	  MatrixXd keyZM;
	  vector<int> keyframes;
	  if(inf.read("keyframes",keyframes) && inf.readMatFile("keyZ",keyZM)){
		vector<VectorXd> keyZ;
		EIGEN3EXT::convert(keyZM,keyZ);
		interpolator->setKeyframe(keyZ,keyframes);
	  }
	}
  }

  if(_animation)
	_animation->setTotalFrameNum(totalFrameNum());
  return succ;
}

// Animation
int AniEditDM::totalFrameNum()const{

  if (interpolator != NULL){
	return interpolator->getT();
  }
  return 0;
}

pObjmesh_const AniEditDM::getOutputObjMesh(const int frame)const{

  if (frame >= 0 && interpolator != NULL && totalFrameNum() > 0){
	const VectorXd &vol_u = getVolFullU(frame);
	vol_obj->interpolate(vol_u);
  }
  return vol_obj->getObjMesh();
}

pObjmesh_const AniEditDM::getInputObjMesh(const int frame)const{
  
  if (frame >= 0 && interpolator != NULL && totalFrameNum() > 0){
	const VectorXd &vol_u = getInputU(frame);
	vol_obj->interpolate(vol_u);
  }
  return vol_obj->getObjMesh();
}

pObjmesh_const AniEditDM::getKeyframeObjMesh()const{
  
  pObjmesh_const pobj;
  const int f = currentFrameNum();
  if (interpolator && interpolator->isKeyframe(f)){
	const VectorXd& key_u_full = interpolator->getKeyframe(f);
	if (key_u_full.size() > 0){
	  vol_obj->interpolate(key_u_full);
	  pobj = vol_obj->getObjMesh();
	}
  }
  return pobj;
}

const VectorXd &AniEditDM::getVolFullU(const int frame)const{

  static VectorXd tempt_u0;
  assert_in (frame, 0, totalFrameNum() );
  if (interpolator != NULL && totalFrameNum() > 0){
	return interpolator->getInterpU(frame);
  }else{
	if (tempt_u0.size() != fullDim()){
	  tempt_u0.resize(fullDim());
	  tempt_u0.setZero();
	}
	return tempt_u0;
  }
}

const VectorXd &AniEditDM::getUnwarpedU()const{
  
  static VectorXd tempt_u0;
  if (interpolator != NULL && totalFrameNum() > 0){
	return interpolator->getUnwarpedU(currentFrameNum());
  }else{
	if (tempt_u0.size() != fullDim()){
	  tempt_u0.resize(fullDim());
	  tempt_u0.setZero();
	}
	return tempt_u0;
  }
}

const VectorXd &AniEditDM::getUforConstraint()const{
  
  static VectorXd tempt_u0;
  if (interpolator != NULL && totalFrameNum() > 0){
	return interpolator->getUforConstraint(currentFrameNum());
  }else{
	if (tempt_u0.size() != fullDim()){
	  tempt_u0.resize(fullDim());
	  tempt_u0.setZero();
	}
	return tempt_u0;
  }
}

const VectorXd &AniEditDM::getInputU(const int frame)const{
  
  static VectorXd tempt_u0;
  assert_in (frame, 0, totalFrameNum() );
  if (interpolator != NULL && totalFrameNum() > 0){
	return interpolator->getInputU(frame);
  }else{
	if (tempt_u0.size() != fullDim()){
	  tempt_u0.resize(fullDim());
	  tempt_u0.setZero();
	}
	return tempt_u0;
  }
}

vector<set<int> > AniEditDM::getConNodes()const{

  const int frame_id = currentFrameNum();
  pConNodesOfFrame_const con = con_nodes_for_edit.getConNodeGroup(frame_id);
  if (con != NULL){
	return con->getConNodesSet();
  }else{
	vector<set<int> > empty;
	return empty;
  }
}

// constraint nodes
void AniEditDM::addConNodes(const vector<int> &group){
  
  if( group.size()>0 && this->editable()){
	
  	const int frame_id = currentFrameNum();
  	con_nodes_for_edit.addConNodeGroup(group, frame_id);
	pConNodesOfFrame_const pc = con_nodes_for_edit.getConNodeGroup(frame_id);
	assert (pc != NULL);
	const VectorXd bary_rest = barycenOfRestShape(pc->getConNodesSet());
	const VectorXd bary_uc = barycenOfVolU(pc->getConNodesSet());
	con_nodes_for_edit.setConPos(bary_uc,bary_rest,frame_id);
	if (interpolator != NULL){
	  const vector<set<int> > &groups = pc->getConNodesSet();
	  interpolator->setConGroups(currentFrameNum() ,groups, bary_uc);
	}
  }
}

void AniEditDM::rmConNodes(const vector<int> &group){
  
  if( group.size()>0 ){

  	const int frame_id = currentFrameNum();
  	con_nodes_for_edit.rmConNodeGroup(group, frame_id);
	pConNodesOfFrame_const pc = con_nodes_for_edit.getConNodeGroup(frame_id);
	if (pc != NULL){
	  const VectorXd bary_rest = barycenOfRestShape(pc->getConNodesSet());
	  const VectorXd bary_uc = barycenOfVolU(pc->getConNodesSet());
	  con_nodes_for_edit.setConPos(bary_uc,bary_rest,frame_id);
	  if (interpolator != NULL){
		const vector<set<int> > &groups = pc->getConNodesSet();
		interpolator->setConGroups(currentFrameNum() ,groups, bary_uc);
	  }
	}
  }
}

void AniEditDM::updateConPos(const VectorXd &uc){

  con_nodes_for_edit.updateConNodeGroup(uc,currentFrameNum());
  if (interpolator != NULL && totalFrameNum() > 0){
	recordDragOp(currentFrameNum(), uc);
	interpolator->setUc(currentFrameNum(), uc);
	this->interpolate();
  }
}

void AniEditDM::removeAllPosCon(){

  con_nodes_for_warping.clear();
  con_nodes_for_edit.clear();
  if (interpolator != NULL){
	interpolator->removeAllPosCon();
  }
}

// interpolate
bool AniEditDM::interpolate(){

  bool succ = false;
  if (interpolator != NULL){
	
	UTILITY::Timer timer;
	timer.start();
	succ = interpolator->interpolate();
	timer.stop("total interpolate time: ");
  }
  if (!succ){
	ERROR_LOG("interpolation failed");
  }
  return succ;
}

// print information
void AniEditDM::print()const{

  if(vol_obj){
	INFO_LOG( "obj vertex number " << vol_obj->getObjMesh()->getVertsNum());
	INFO_LOG( "vol vertex number " << vol_obj->getTetMesh()->nodes().size());
	INFO_LOG( "vol tet number " << vol_obj->getTetMesh()->tets().size());
  }
  INFO_LOG( "reduce dimension " << this->reducedDim());
  INFO_LOG( "total frame number " << this->totalFrameNum());
}

// file io
bool AniEditDM::saveParitalCon(const string filename)const{

  return con_nodes_for_edit.save(filename);
}

bool AniEditDM::loadParitalCon(const string filename){
  
  const bool succ = con_nodes_for_edit.load(filename);
  if (interpolator != NULL && succ){
	interpolator->removeAllPosCon();
	const set<pConNodesOfFrame>&con_groups=con_nodes_for_edit.getConNodeGroups();
	BOOST_FOREACH(pConNodesOfFrame pc, con_groups){
	  if ( pc != NULL && !pc->isEmpty() ) {
		const int frame_id = pc->getFrameId();
		if(_animation)
		  _animation->setCurrentFrame(frame_id);
		const VectorXd bary_uc = pc->getBarycenterUc();
		const vector<set<int> > &groups = pc->getConNodesSet();
		interpolator->setConGroups(frame_id,groups, bary_uc);
	  }
	}
  }
  return succ;
}

string AniEditDM::getZeroStr(const int frame, const int T)const{
  
  if (frame < 10){
	return string("000");
  }else if (frame < 100){
	return string("00");
  }else if (frame < 1000){
	return string("0");
  }else{
	// @bug
	ERROR_LOG("the frame number is too large: " << frame);
	return string("_0");
  }
}

bool AniEditDM::saveOutputMeshes(const string filename)const{

  bool succ = false;
  const int T = totalFrameNum();
  for (int f = 0; f < T; ++f){

  	const string f_id = boost::lexical_cast<string>(f);
    const string frame = filename +getZeroStr(f,T)+ f_id + string(".obj");
  	cout << frame << endl;
  	succ = getOutputObjMesh(f)->write(frame);
  	if(!succ){
  	  ERROR_LOG("failed to save frame " << f << " to " << frame);
  	  break;
  	}
  }
  return succ;
}

bool AniEditDM::saveInputMeshes(const string filename)const{

  bool succ = false;
  const int T = totalFrameNum();
  for (int f = 0; f < T; ++f){

  	const string f_id = boost::lexical_cast<string>(f);
    const string frame = filename +getZeroStr(f,T)+ f_id + string(".obj");
	INFO_LOG(frame);
	succ = getInputObjMesh(f)->write(frame);
  	if(!succ){
  	  ERROR_LOG("failed to save frame " << f << " to " << frame);
  	  break;
  	}
  }
  return succ;
}

bool AniEditDM::saveCurrentOutputMesh(const string filename)const{
  return getOutputObjMesh(currentFrameNum())->write(filename);
}

bool AniEditDM::saveCurrentInputMesh(const string filename)const{
  return getInputObjMesh(currentFrameNum())->write(filename);
}

bool AniEditDM::saveCurrentOutputVolMesh(const string filename){
  
  pTetMesh vol_mesh = getVolMesh();
  bool succ = false;
  if (vol_mesh) {
	const VectorXd &u = getVolFullU(currentFrameNum());
	if(u.size() > 0){
	  const VectorXd _u = -1.0f*u;
	  vol_mesh->applyDeformation(u);
	  succ = vol_mesh->write(filename);
	  vol_mesh->applyDeformation(_u);
	}
  }
  return succ;
}

bool AniEditDM::saveVolFullU(const string filename)const{
  
  const int T = totalFrameNum();
  const int n = fullDim();
  MatrixXd U(n,T);
  for (int i = 0; i < T; ++i){
    U.col(i) = getVolFullU(i);
  }
  return EIGEN3EXT::write(filename,U);
}

bool AniEditDM::saveCurrentVolFullU(const string filename)const{
  
  const VectorXd &u = getVolFullU();
  return EIGEN3EXT::write(filename,u);
}

bool AniEditDM::loadCurrentReducedEdits(const string filename){

  bool succ = false;
  const int f = currentFrameNum();
  if (interpolator && f >= 0 && f < interpolator->getT()){
	VectorXd z;
	succ = EIGEN3EXT::load(filename,z);
	if(succ){
	  assert_eq(z.size(),interpolator->reducedDim());
	  interpolator->setReducedEdits(f,z);
	}
  }
  return succ;
}

bool AniEditDM::saveCurrentReducedEdits(const string filename)const{
  
  bool succ = false;
  const int f = currentFrameNum();
  if (interpolator && f >= 0 && f < interpolator->getT()){
	const VectorXd& z = interpolator->getReducedEdits(f);
	succ = EIGEN3EXT::write(filename,z);
  }
  return succ;
}

bool AniEditDM::loadAllReducedEdits(const string filename)const{
  
  vector<VectorXd> z;
  bool succ = EIGEN3EXT::load(filename,z);
  if(succ){
	assert_eq((int)z.size(),interpolator->getT());
	if(z.size() > 0){
	  assert_eq(z[0].size(),interpolator->reducedDim());
	}
	interpolator->setReducedEdits(z);
  }
  return succ;
}

bool AniEditDM::saveAllReducedEdits(const string filename)const{

  assert(interpolator);
  assert_eq(interpolator->getT(), totalFrameNum());
  MatrixXd allZ(reducedDim(),totalFrameNum());
  for (int f = 0; f < totalFrameNum(); ++f){
	allZ.col(f) = interpolator->getReducedEdits(f);
  }
  return EIGEN3EXT::load(filename,allZ);
}

// return the barycenters of a set of nodes of the rest volume shape.
VectorXd AniEditDM::barycenOfRestShape(const vector<set<int> >&g)const{

  VectorXd rest_barycenters(g.size()*3);
  rest_barycenters.setZero();
  pTetMesh_const tet = getVolMesh();
  assert(tet);
  const VVec3d &nodes = tet->nodes();
  for (size_t i = 0; i < g.size(); ++i){
	const set<int> &one_group = g[i];
	set<int>::const_iterator it = one_group.begin();
	for (; it != one_group.end(); ++it)
	  rest_barycenters.segment(i*3,3) += nodes[*it];
	if (one_group.size() > 0)
	  rest_barycenters.segment(i*3,3) /= one_group.size();
  }
  return rest_barycenters;
}

// return barycenters of vol_u.
VectorXd AniEditDM::barycenOfVolU(const vector<set<int> >&g)const{
  
  VectorXd uc(g.size()*3);
  uc.setZero();
  if (currentFrameNum() >=0){
	const VectorXd &vol_u = this->getVolFullU();
	for (size_t i = 0; i < g.size(); ++i){
	
	  BOOST_FOREACH(const int v_i, g[i]){
		assert_in (v_i,0,vol_u.size());
		uc[i*3+0] = vol_u[v_i*3+0];
		uc[i*3+1] = vol_u[v_i*3+1];
		uc[i*3+2] = vol_u[v_i*3+2];
	  }
	  if (g[i].size() > 0){

		uc[i*3+0] /= g[i].size();
		uc[i*3+1] /= g[i].size();
		uc[i*3+2] /= g[i].size();
	  }
	}
  }
  return uc;
}

bool AniEditDM::editable()const{
  
  if (interpolator != NULL){
	return interpolator->editable(currentFrameNum());
  }
  return false;
}

bool AniEditDM::loadConPath(const string filename){

  TRACE_FUN();

  // read in the constrained nodes and frames' ids.
  vector<int> con_frame_ids;
  vector<set<int> > con_nodes(1);
  vector<VectorXd> target_pos;
  ifstream in_file;
  in_file.open (filename.c_str());
  if (!in_file.is_open()){
	ERROR_LOG("failed to open file: "<<filename);
	return false;
  }

  // load constrained nodes
  int con_node_size = 0;
  in_file >> con_node_size;
  assert_gt(con_node_size,0);
  con_nodes.resize(con_node_size);
  for (int i = 0; i < con_node_size; ++i){
	int con_node_id = 0;
    in_file >> con_node_id;
	con_nodes[0].insert(con_node_id);;
  }

  // load constrained frames and target positions
  int con_frame_num = 0;
  in_file >> con_frame_num;
  assert_gt(con_node_size,0);
  con_frame_ids.resize(con_frame_num);
  target_pos.resize(con_frame_num);
  for (int i = 0; i < con_frame_num; ++i){
	target_pos[i].resize(3);
    in_file >> con_frame_ids[i];
    in_file >> target_pos[i][0];
    in_file >> target_pos[i][1];
    in_file >> target_pos[i][2];
  }

  // set the constrained nodes
  assert_eq(con_frame_ids.size(), target_pos.size());
  con_nodes_for_edit.clear();
  const VectorXd bary_rest = this->barycenOfRestShape(con_nodes);
  for (size_t i=0; i < con_frame_ids.size(); ++i){
	pConNodesOfFrame con_frame = pConNodesOfFrame(new ConNodesOfFrame());
	const VectorXd bary_uc = target_pos[i]-bary_rest;
	con_frame->setConNodes(con_nodes,bary_uc,bary_rest,con_frame_ids[i]);
	con_nodes_for_edit.addConNodeGroup(con_frame);
	if (interpolator != NULL)
	  interpolator->setConGroups(con_frame_ids[i] ,con_nodes, bary_uc);
  }
  return true;
}

int AniEditDM::reducedDim()const{

  if (interpolator){
	return interpolator->reducedDim();
  }else{
	return 0;
  }
}

int AniEditDM::fullDim()const{
  if (vol_obj){
	return vol_obj->getTetMesh()->nodes().size()*3;
  }
  return 0;
}

bool AniEditDM::playRecodDrag(){
  
  bool succ = true;
  this->interpolate();
  for (int i = 0; i < drag_record.totalRecord(); ++i){

	const OneDragRecord &one_drag = drag_record.getRecord(i);
	interpolator->setUc(one_drag.getFrameNum(), one_drag.getUc());
	if ( !this->interpolate() ){
	  succ = false;
	  break;
	}
  }
  return succ;
}

void AniEditDM::computeConNodeTrajectory(){
  
  pTetMesh_const vol_mesh = this->getVolMesh();
  const vector<set<int> > groups = this->getConNodes();
  size_t total_con_node = 0;
  for (size_t i = 0; i < groups.size(); ++i){
	total_con_node += groups[i].size();
  }
  con_node_traj.resize(total_con_node);
  for (size_t i = 0; i < con_node_traj.size(); ++i){
    con_node_traj[i].resize(3*totalFrameNum());
	con_node_traj[i].setZero();
  }
  for (int f = 0; f < totalFrameNum(); ++f){
	const VectorXd &vol_u = this->getVolFullU(f);
	int count = 0;
	for (size_t g = 0; g < groups.size(); ++g){
	  BOOST_FOREACH(int i, groups[g]){
		const Vector3d &rest_v = vol_mesh->nodes()[i];
		con_node_traj[count][f*3+0] = vol_u[i*3+0] + rest_v[0];
		con_node_traj[count][f*3+1] = vol_u[i*3+1] + rest_v[1];
		con_node_traj[count][f*3+2] = vol_u[i*3+2] + rest_v[2];
		count ++;
	  }
	}
  }
}
