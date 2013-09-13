#include <VolBBox.h>
#include "VolObjMeshCtrl.h"
using namespace LSW_RENDERMESH;
using namespace LSW_BASE_UI;

VolObjMeshCtrl::VolObjMeshCtrl(QMainWindow *main_win, pVolObjMesh volobjmesh):
  volobjmesh(volobjmesh){
  
  this->file_dialog = pFileDialog(new FileDialog(main_win));
}

bool VolObjMeshCtrl::initialize(const string filename){
  
  if (volobjmesh == NULL)
	volobjmesh = pVolObjMesh(new VolObjMesh());
  const bool succ = volobjmesh->initialize(filename);
  sendMsgResetScene();
  return succ;
}

void VolObjMeshCtrl::loadObjMesh(){

  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(volobjmesh != NULL && volobjmesh->loadObjMesh(filename));
	sendMsgResetScene();
  }
}

void VolObjMeshCtrl::loadVolSurface(){

  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(volobjmesh != NULL && volobjmesh->loadVolSurface(filename));
	sendMsgResetScene();
  }
}

void VolObjMeshCtrl::loadVolMesh(){

  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(volobjmesh != NULL && volobjmesh->loadVolMesh(filename));
	sendMsgResetScene();
  }
}

void VolObjMeshCtrl::loadInterpWeights(){
	  
  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(volobjmesh != NULL && volobjmesh->loadInterpWeights(filename));
  }
}

void VolObjMeshCtrl::saveInterpWeights()const{

  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(volobjmesh != NULL && volobjmesh->saveInterpWeights(filename));
  }
}

void VolObjMeshCtrl::precomputeInterpWeights(){

  bool succ = false;
  if(volobjmesh != NULL){
	succ = volobjmesh->precomputeInterpWeights();
  }
  file_dialog->warning(succ,"failed to perferm mesh interpolation!");
}

void VolObjMeshCtrl::sendMsgResetScene()const{
  
  if(this->getVolObjMesh() != NULL && this->getVolMesh() != NULL){

	double min[3],max[3];
	VolBBox::getInstance()->getBBbox(this->getVolMesh(),min,max);
	emit resetSceneMsg(min[0],min[1],min[2],max[0],max[1],max[2]);
  }
}

void VolObjMeshCtrl::togglePhoneShading(){
  
  if (volobjmesh && volobjmesh->getObjMesh()){
	const bool use = volobjmesh->getObjMesh()->usePhoneShading();
	volobjmesh->getObjMesh()->setPhoneShading(!use);
	emit update();
  }
}
