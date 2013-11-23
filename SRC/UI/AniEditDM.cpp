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
  _animation(animation), vol_obj(vol_obj){}

AniEditDM::AniEditDM(pTetMeshEmbeding vol_obj, pAniDataModel animation, pBaseInterpolator interpolator):
  _animation(animation),interpolator(interpolator), vol_obj(vol_obj){}

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
	  if ( inf.readFilePath("partial_constraints",partial_con) )
		succ &= loadParitalCon(partial_con);
	  
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

  pPartialConstraints_const con = _partialCon.getPartialCon(currentFrameNum());
  if (con) return con->getConNodesSet();
  return vector<set<int> >();
}

// constraint nodes
void AniEditDM::addConNodes(const vector<int> &group){
  if (group.size()>0 && editable() && interpolator){
	_partialCon.addConNodes(group,currentFrameNum());
	resetPartialCon(currentFrameNum());
  }
}

void AniEditDM::rmConNodes(const vector<int> &group){
  if (group.size()>0 && editable() && interpolator){
	_partialCon.rmConNodes(group,currentFrameNum());
	resetPartialCon(currentFrameNum());
  }
}

void AniEditDM::updateConPos(const Vector3d &uc,const int group_id){
  pPartialConstraints par = getPartialCon();
  if (par)	par->moveOneGroup(uc,group_id);
}

void AniEditDM::updateConPos(const VectorXd &uc){
  if(interpolator != NULL && totalFrameNum() > 0){
	if(_partialCon.updatePc(uc,currentFrameNum())){
	  interpolator->setUc(currentFrameNum(), uc);
	  this->interpolate();
	}
  }
}

void AniEditDM::removeAllPosCon(){
  _partialCon.clear();
  if(interpolator){
	interpolator->removeAllPosCon();
  }
}

void AniEditDM::resetPartialCon(const int frameid){
  
  pPartialConstraints_const par = _partialCon.getPartialCon(frameid);
  if (par){

	const VectorXd &u = getUforConstraint();
	Matrix<double,3,-1> pc(3,par->numConNodes());
	const vector<set<int> > &vs = par->getConNodesSet();
	BOOST_FOREACH(const set<int>& s, vs){
	  BOOST_FOREACH(const int i, s){
		assert_in(i,0,pc.cols()-1);
		assert_in(i*3,0,u.size()-3);
		pc.col(i) = u.segment<3>(i*3);
	  }
	}
	_partialCon.updatePc(pc,frameid);
	interpolator->setConGroups(par->getFrameId(),par->getConNodesSet(),par->getPc());
  }else{
	interpolator->removePartialCon(frameid);
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

  return _partialCon.save(filename);
}

bool AniEditDM::loadParitalCon(const string filename){
  
  const bool succ = _partialCon.load(filename);
  if (interpolator != NULL && succ){
	interpolator->removeAllPosCon();
	const set<pPartialConstraints> &con_groups=_partialCon.getPartialConSet();
	BOOST_FOREACH(pPartialConstraints par, con_groups)
	  interpolator->setConGroups(par->getFrameId(),par->getConNodesSet(),par->getPc());
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

bool AniEditDM::saveOutputVolMeshesVTK(const string filename)const{

  pTetMesh_const vol_mesh = getVolMesh();
  bool succ = false;
  if (vol_mesh){
	succ = true;
	for (int i = 0; i < totalFrameNum() && succ; ++i)
	  succ = vol_mesh->writeVTK(filename+TOSTR(i)+".vtk",getVolFullU(i));
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
