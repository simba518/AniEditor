#include <fstream>
#include <Timer.h>
#include <boost/foreach.hpp>
#include <MatrixIO.h>
#include <AuxTools.h>
#include <JsonFilePaser.h>
#include <InterpolatorFactory.h>
#include <MatrixTools.h>
#include "AniEditDM.h"
using namespace std;
using namespace UTILITY;
using namespace ANI_EDIT_UI;

AniEditDM::AniEditDM(pTetMeshEmbeding vol_obj, pAniDataModel animation):
  _animation(animation), vol_obj(vol_obj){
  enableInterpolate = true;
}

AniEditDM::AniEditDM(pTetMeshEmbeding vol_obj, pAniDataModel animation, pBaseInterpolator interpolator):
  _animation(animation),interpolator(interpolator), vol_obj(vol_obj){
  enableInterpolate = true;
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

	  string partial_con, drag_trajectory;
	  if ( inf.readFilePath("partial_con",partial_con) )
		succ &= loadParitalCon(partial_con);

	  if (inf.readFilePath("path_con", partial_con))
		succ &= loadConPath(partial_con);
	  
	  MatrixXd keyZM;
	  vector<int> keyframes;
	  if(inf.read("keyframes",keyframes) && inf.readMatFile("keyZ",keyZM)){
		vector<VectorXd> keyZ;
		EIGEN3EXT::convert(keyZM,keyZ);
		interpolator->setKeyframe(keyZ,keyframes);
	  }
	  string addani;
	  if(inf.readFilePath("additional_ani",addani)){
		loadAnimation(addani);
	  }

	  string sf;
	  if (inf.readFilePath("scene", sf,true)){
		if (!scene) scene = pObjmesh(new Objmesh());
		scene->load(sf);
	  }
	  if(!inf.readVecFile("translation", scene_translation, UTILITY::TEXT)){
		scene_translation = vector<double>(totalFrameNum()*3,0);
	  }
	  ERROR_LOG_COND("the number of the translation vector is not correct: "<<scene_translation.size(),scene_translation.size()==totalFrameNum()*3);
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

pObjmesh_const AniEditDM::getOutputObjMesh(const int frame, const bool translate)const{

  if (frame >= 0 && interpolator != NULL && totalFrameNum() > 0){
	const VectorXd &vol_u = getVolFullU(frame, translate);
	vol_obj->interpolate(vol_u);
  }
  return vol_obj->getObjMesh();
}

pObjmesh_const AniEditDM::getInputObjMesh(const int frame, const bool translate)const{
  
  if (frame >= 0 && interpolator != NULL && totalFrameNum() > 0){
	const VectorXd &vol_u = getInputU(frame, translate);
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

const VectorXd &AniEditDM::getVolFullU(const int frame, const bool translate)const{

  static VectorXd tempt_u0;
  assert_in (frame, 0, totalFrameNum() );
  if (interpolator != NULL && totalFrameNum() > 0){
	if (!translate){
	  return interpolator->getInterpU(frame);
	}else{
	  tempt_u0 = interpolator->getInterpU(frame);
	  this->translate(tempt_u0, frame);
	  return tempt_u0;
	}
  }else{
	if (tempt_u0.size() != fullDim()){
	  tempt_u0.resize(fullDim());
	  tempt_u0.setZero();
	}
	return tempt_u0;
  }
}

const VectorXd &AniEditDM::getUnwarpedU(const int frame, const bool translate)const{
  
  static VectorXd tempt_u0;
  if (interpolator != NULL && totalFrameNum() > 0){
	assert_in(frame,0, totalFrameNum()-1);
	if (!translate){
	  return interpolator->getUnwarpedU(frame);
	}else{
	  tempt_u0 = interpolator->getUnwarpedU(frame);
	  this->translate(tempt_u0, frame);
	  return tempt_u0;
	}
  }else{
	if (tempt_u0.size() != fullDim()){
	  tempt_u0.resize(fullDim());
	  tempt_u0.setZero();
	}
	return tempt_u0;
  }
}

const VectorXd &AniEditDM::getUforConstraint(const int frame)const{
  
  static VectorXd tempt_u0;
  if (interpolator != NULL && totalFrameNum() > 0){
	assert_in(frame,0, totalFrameNum()-1);
	return interpolator->getUforConstraint(frame);
  }else{
	if (tempt_u0.size() != fullDim()){
	  tempt_u0.resize(fullDim());
	  tempt_u0.setZero();
	}
	return tempt_u0;
  }
}

const VectorXd &AniEditDM::getInputU(const int frame, const bool translate)const{
  
  static VectorXd tempt_u0;
  assert_in (frame, 0, totalFrameNum() );
  if (interpolator != NULL && totalFrameNum() > 0){
	if (!translate){
	  return interpolator->getInputU(frame);
	}else{
	  tempt_u0 = interpolator->getInputU(frame);
	  this->translate(tempt_u0, frame);
	  return tempt_u0;
	}
  }else{
	if (tempt_u0.size() != fullDim()){
	  tempt_u0.resize(fullDim());
	  tempt_u0.setZero();
	}
	return tempt_u0;
  }
}

vector<set<int> > AniEditDM::getConNodes()const{

  pPartialConstraints_const con = _partialCon.getPartialCon(currentFrameNum());
  if (con) return con->getConNodesSet();
  return vector<set<int> >();
}

// constraint nodes
const Matrix<double,3,-1> AniEditDM::getXc(const int group)const{
  
  pPartialConstraints_const _partialCon = getPartialCon();
  pPartialConstraints_const _partialCon_x0 = getPartialConX0();
  assert(_partialCon);
  assert(_partialCon_x0);
  assert_in(group,0,_partialCon->numGroup()-1);
  const Matrix<double,3,-1> xc = _partialCon_x0->getPc(group);
  const Matrix<double,3,-1> uc = _partialCon->getPc(group);
  assert_eq(xc.size(), uc.size());
  return xc + uc;
}

void AniEditDM::updateXc(const Matrix<double,3,-1> &xc,const int group_id){

  pPartialConstraints_const _partialCon_x0 = getPartialConX0();
  assert(_partialCon_x0);
  assert_in(group_id,0,_partialCon_x0->numGroup()-1);
  const Matrix<double,3,-1> xc_0 = _partialCon_x0->getPc(group_id);
  assert_eq(xc_0.size(), xc.size());
  const Matrix<double,3,-1> uc = xc - xc_0;
  updateConPos(uc, group_id);
}

void AniEditDM::addPartialCon(const PartialConstraintsSet&parcons){

  _partialCon.addPartialCon(parcons);
  const set<pPartialConstraints> &cons = _partialCon.getPartialConSet();
  BOOST_FOREACH(pPartialConstraints ele, cons){
	assert_eq(ele->getPc().size(),ele->numConNodes()*3);
	interpolator->setConGroups(ele->getFrameId(),ele->getConNodesSet(),ele->getPc());
  }
  // this->interpolate();
}

void AniEditDM::addConNodes(const vector<int> &group){
  if (group.size()>0 && editable() && interpolator){
	_partialCon.addConNodes(group,currentFrameNum());
	_partialCon_x0.addConNodes(group,currentFrameNum());
	resetPartialCon(currentFrameNum());
  }
}

void AniEditDM::rmConNodes(const vector<int> &group){
  if (group.size()>0 && editable() && interpolator){
	_partialCon.rmConNodes(group,currentFrameNum());
	_partialCon_x0.rmConNodes(group,currentFrameNum());
	resetPartialCon(currentFrameNum());
  }
}

void AniEditDM::updateConPos(const Matrix<double,3,-1> &uc,const int group_id){

  pPartialConstraints par = getPartialCon();
  if (par){
	par->updatePc(uc,group_id);
	if(interpolator != NULL && totalFrameNum() > 0){
	  interpolator->setUc(currentFrameNum(), par->getPc());
	  this->interpolate();
	}
  }
}

void AniEditDM::updateConPos(const Matrix<double,3,-1> &uc){
  if(interpolator != NULL && totalFrameNum() > 0){
	if(_partialCon.updatePc(uc,currentFrameNum())){
	  interpolator->setUc(currentFrameNum(), uc);
	  this->interpolate();
	}
  }
}

void AniEditDM::removeAllPosCon(){
  _partialCon.clear();
  _partialCon_x0.clear();
  if(interpolator){
	interpolator->removeAllPosCon();
  }
}

void AniEditDM::resetPartialCon(const int frameid, const bool interp){
  
  pPartialConstraints_const par = _partialCon.getPartialCon(frameid);
  if (par){

	const VectorXd &u = getUforConstraint(frameid);
	Matrix<double,3,-1> pc;
	getSubUc(par->getConNodesSet(),u,pc);
	_partialCon.updatePc(pc,frameid);
	interpolator->setConGroups(par->getFrameId(),par->getConNodesSet(),par->getPc());

	const pTetMesh_const tetmesh = getVolMesh();
	if (tetmesh && rest_shape.size() != tetmesh->nodes().size()*3){
	  tetmesh->nodes(rest_shape);
	}
	assert_eq(getInputU(frameid).size(), rest_shape.size());
	const VectorXd ref_x = getInputU(frameid)+rest_shape;
	getSubUc(par->getConNodesSet(),ref_x,pc);
	_partialCon_x0.updatePc(pc,frameid);
	
  }else{
	interpolator->removePartialCon(frameid);
  }
  if (interp) this->interpolate();
}

void AniEditDM::getSubUc(const vector<set<int> > &groups,const VectorXd &full_u,Matrix<double,3,-1> &sub_u)const{

  int nodes = 0;
  BOOST_FOREACH(const set<int>& s, groups)
	nodes += s.size();

  sub_u.resize(3,nodes);
  int index = 0;
  BOOST_FOREACH(const set<int>& s, groups){
	BOOST_FOREACH(const int i, s){
	  assert_in(i*3,0,full_u.size()-3);
	  sub_u.col(index) = full_u.segment<3>(i*3);
	  index++;
	}
  }
}

const Matrix<double,3,-1> AniEditDM::getUc(const int group)const{

  pPartialConstraints_const par = _partialCon.getPartialCon(currentFrameNum());
  assert(par);
  return par->getPc(group);
}

bool AniEditDM::loadConPath(const string filename){
  
  bool succ = false;
  JsonFilePaser jsonf;

  if (jsonf.open(filename)){

	vector<int> frame;
	vector<double> pos;
	int con_node_id = 0;
	if(jsonf.read("con_frames", frame) && 
	   jsonf.read("con_pos", pos) &&
	   jsonf.read("con_node_id",con_node_id)){

	  assert_eq(frame.size()*3,pos.size());
	  assert_in(con_node_id,0,vol_obj->getTetMesh()->nodes().size());

	  vector<int> cons;
	  cons.push_back(con_node_id);
	  const Vector3d rc = vol_obj->getTetMesh()->nodes()[con_node_id];

	  Vector3d pc;
	  removeAllPosCon();
	  for (int i = 0; i < frame.size(); ++i){
		pc[0] = pos[i*3+0]-rc[0];
		pc[1] = pos[i*3+1]-rc[1];
		pc[2] = pos[i*3+2]-rc[2];
		_partialCon.addConNodes(cons,frame[i]);
		_partialCon.updatePc(pc,frame[i]);
		pPartialConstraints_const par = _partialCon.getPartialCon(frame[i]);
		assert(par);
		interpolator->setConGroups(par->getFrameId(),par->getConNodesSet(),par->getPc());

		_partialCon_x0.addConNodes(cons,frame[i]);
		_partialCon_x0.updatePc(rc,frame[i]);
	  }
	  succ = true;
	}
  }

  ERROR_LOG_COND("failed to load data from path file: " << filename, succ);
  return succ;
}

// interpolate
bool AniEditDM::interpolate(){

  if (!enableInterpolate)
	return true;

  bool succ = false;
  if (interpolator != NULL){
	
	// UTILITY::Timer timer;
	// timer.start();
	succ = interpolator->interpolate();
	// timer.stop("total interpolate time: ");
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
  if(interpolator){
	interpolator->printSolveInf();
  }
  printConPositions();
}

// file io
bool AniEditDM::saveParitalCon(const string filename)const{
  return _partialCon.write(filename);
}

bool AniEditDM::loadParitalCon(const string filename){
  
  const bool succ = _partialCon.load(filename);
  if (interpolator != NULL && succ){
	interpolator->removeAllPosCon();
	const set<pPartialConstraints> &con_groups=_partialCon.getPartialConSet();
	BOOST_FOREACH(pPartialConstraints par, con_groups)
	  interpolator->setConGroups(par->getFrameId(),par->getConNodesSet(),par->getPc());
  }

  if(succ){

	_partialCon_x0.load(filename);
	set<pPartialConstraints> &con_groups=_partialCon_x0.getPartialConSet();
	Matrix<double,3,-1> pc;
	BOOST_FOREACH(pPartialConstraints par, con_groups){

	  const int frameid = par->getFrameId();
	  const pTetMesh_const tetmesh = getVolMesh();
	  if (tetmesh && rest_shape.size() != tetmesh->nodes().size()*3)
		tetmesh->nodes(rest_shape);
	  assert_eq(getInputU(frameid).size(), rest_shape.size());
	  const VectorXd ref_x = getInputU(frameid)+rest_shape;
	  getSubUc(par->getConNodesSet(),ref_x,pc);
	  _partialCon_x0.updatePc(pc,frameid);
	}
  }
  return succ;
}

void AniEditDM::printConPositions()const{
  
  const VVec3d& rest_u = getVolMesh()->nodes();
  const set<pPartialConstraints>& parcon = _partialCon.getPartialConSet();

  BOOST_FOREACH(pPartialConstraints_const con, parcon){

	if (con){
	  const vector<set<int> > &group = con->getConNodesSet();
	  const Matrix<double,3,-1> &con_u = con->getPc();
	  int p = 0;
	  cout << "frame: " << con->getFrameId() << endl;
	  for (int i = 0; i < group.size(); ++i){

		BOOST_FOREACH(const int &node_i, group[i]){
		  assert_in(p,0,con_u.cols()-1);
		  assert_in(node_i,0,rest_u.size()-1);
		  const Vector3d v = con_u.col(p)+rest_u[node_i];
		  cout<< v[0]<< ", " << v[1]<< ", " << v[2]<< ", " << endl;
		  p ++;
		}
	  }
	}
  }
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
  	succ = getOutputObjMesh(f,true)->write(frame);
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
	succ = getInputObjMesh(f,true)->write(frame);
  	if(!succ){
  	  ERROR_LOG("failed to save frame " << f << " to " << frame);
  	  break;
  	}
  }
  return succ;
}

bool AniEditDM::saveOutputVolMeshesVTK(const string filename)const{

  pTetMesh_const vol_mesh = getVolMesh();
  bool succ = false;
  if (vol_mesh){
	succ = true;
	for (int i = 0; i < totalFrameNum() && succ; ++i)
	  succ = vol_mesh->writeVTK(filename+TOSTR(i)+".vtk",getVolFullU(i,true));
  }
  return succ;
}

bool AniEditDM::saveOutputObjMeshesVTK(const string filename)const{
  
  pObjmesh_const obj_mesh = getOutputObjMesh();
  assert(obj_mesh);

  bool succ = false;
  for (int i = 0; i < totalFrameNum(); ++i){
	succ = getOutputObjMesh(i,true)->writeVTK(filename+TOSTR(i)+".vtk");
	if (!succ) break;
  }
  return succ;
}

bool AniEditDM::saveAdditionalAniObjMeshesVTK(const string filename)const{
  
  bool succ = false;
  for (int i = 0; i < totalFrameNum(); ++i){
	if (hasAdditionalAni()){
	  const VectorXd &vol_u = additionalVolU[i];
	  vol_obj->interpolate(vol_u);
	  succ = vol_obj->getObjMesh()->writeVTK(filename+TOSTR(i)+".vtk");
	  if (!succ) break;
	}
  }
  return succ;
}

bool AniEditDM::saveAdditionalAniObjMeshes(const string filename)const{

  bool succ = false;
  const int T = totalFrameNum();
  for (int i = 0; i < totalFrameNum(); ++i){
	if (hasAdditionalAni()){
	  const string f_id = boost::lexical_cast<string>(i);
	  const string frame = filename +getZeroStr(i,T)+ f_id + string(".obj");
	  cout << frame << endl;
	  const VectorXd &vol_u = additionalVolU[i];
	  vol_obj->interpolate(vol_u);
	  succ = vol_obj->getObjMesh()->write(frame);
	  if (!succ) break;
	}
  }
  return succ;
}

bool AniEditDM::saveCurrentOutputMesh(const string filename)const{
  return getOutputObjMesh(currentFrameNum(),true)->write(filename);
}

bool AniEditDM::saveCurrentInputMesh(const string filename)const{
  return getInputObjMesh(currentFrameNum(),true)->write(filename);
}

bool AniEditDM::saveCurrentOutputVolMesh(const string filename){
  
  pTetMesh vol_mesh = getVolMesh();
  bool succ = false;
  if (vol_mesh) {
	const VectorXd &u = getVolFullU(currentFrameNum(),true);
	if(u.size() > 0){
	  const VectorXd _u = -1.0f*u;
	  vol_mesh->applyDeformation(u);
	  succ = vol_mesh->write(filename);
	  vol_mesh->applyDeformation(_u);
	}
  }
  return succ;
}

bool AniEditDM::saveSceneSequence(const string filename, const bool VTK)const{

  if (!scene){
	ERROR_LOG("there is no scene file.");
	return false;
  }
  bool succ = false;
  const int T = totalFrameNum();
  Objmesh t_scene = *scene;
  for (int f = 0; f < T; ++f){

  	const string f_id = boost::lexical_cast<string>(f);
	const string append = VTK? string(".vtk"):string(".obj");
    const string frame = filename +getZeroStr(f,T)+ f_id + append;
  	cout << frame << endl;
	VectorXd v = scene->getVerts();
	this->translate(v,f);
	t_scene.setVerts(v);
	if (VTK){
	  succ = t_scene.writeVTK(frame);
	}else{
	  succ = t_scene.write(frame);
	}
  	if(!succ){
  	  ERROR_LOG("failed to save frame " << f << " to " << frame);
  	  break;
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

bool AniEditDM::loadAllReducedEdits(const string filename){
  
  vector<VectorXd> z;
  bool succ = EIGEN3EXT::load(filename,z);
  if(succ){
	assert_ge((int)z.size(),interpolator->getT());
	if(z.size() > 0) assert_eq(z[0].size(),interpolator->reducedDim());
	if (z.size() > interpolator->getT()){
	  vector<VectorXd> z_sub;
	  for (int i = 0; i < interpolator->getT(); ++i)
		z_sub.push_back(z[i]);
	  interpolator->setReducedEdits(z_sub);
	}else{
	  interpolator->setReducedEdits(z);
	}
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
  return EIGEN3EXT::write(filename,allZ);
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

int AniEditDM::reducedDim()const{

  if (interpolator){
	return interpolator->reducedDim();
  }else{
	return 0;
  }
}

int AniEditDM::fullDim()const{
  return vol_obj ? vol_obj->getTetMesh()->nodes().size()*3:0;
}

void AniEditDM::getSceneTranslate(double &tx,double &ty,double &tz,const int frame)const{

  if (frame*3+2 >= scene_translation.size()){
	tx = ty = tz = 0.0f;
  }else{
	tx = scene_translation[frame*3+0];
	ty = scene_translation[frame*3+1];
	tz = scene_translation[frame*3+2];
  }
  // tx = frame*0.003; // slower dino skate
}

void AniEditDM::translate(VectorXd &u, const int frame)const{

  assert_eq(u.size() % 3,0);
  double x,y,z;
  getSceneTranslate(x,y,z,frame);
  for (int i = 0; i < u.size()/3; ++i){
    u[i*3+0] += x;
    u[i*3+1] += y;
    u[i*3+2] += z;
  }
}

bool AniEditDM::savePartialConBalls(const string filename)const{
  
  Objmesh ball;
  const string ball_f = "./con_ball.obj";
  if (!ball.load(ball_f)){
	ERROR_LOG("failed to load the conball from "<<ball_f);
	return false;
  }

  bool succ = true;

  const VVec3d& rest_u = getVolMesh()->nodes();
  const set<pPartialConstraints>& parcon = getAllConNodes().getPartialConSet();
  BOOST_FOREACH(pPartialConstraints_const con, parcon){

	if (con){

	  Objmesh target_pos_balls;
	  Objmesh input_pos_balls;
	  Objmesh output_pos_balls;

	  const int f = con->getFrameId();
	  const VectorXd &output_u = getVolFullU(f);
	  const VectorXd &input_u = getInputU(f);
	  const vector<set<int> > &group = con->getConNodesSet();
	  const Matrix<double,3,-1> &con_u = con->getPc();
	  Vector3d trans;
	  getSceneTranslate(trans[0], trans[1], trans[2],f);

	  int p = 0;
	  for (int i = 0; i < group.size(); ++i){

		BOOST_FOREACH(const int &node_i, group[i]){

		  assert_in(p,0,con_u.cols()-1);
		  assert_in(node_i,0,rest_u.size()-1);

		  // save target positions
		  ball.moveCenterTo(con_u.col(p)+rest_u[node_i]+trans);
		  target_pos_balls.combine(ball);

		  ball.moveCenterTo(input_u.segment<3>(node_i*3)+rest_u[node_i]+trans);
		  input_pos_balls.combine(ball);

		  ball.moveCenterTo(output_u.segment<3>(node_i*3)+rest_u[node_i]+trans);
		  output_pos_balls.combine(ball);

		  p ++;
		}
	  }
	  succ = target_pos_balls.write(filename+"-target-frame"+TOSTR(f)+".obj");
	  succ = input_pos_balls.write(filename+"-input-frame"+TOSTR(f)+".obj");
	  succ = output_pos_balls.write(filename+"-output-frame"+TOSTR(f)+".obj");
	}
  }

  return succ;
}
