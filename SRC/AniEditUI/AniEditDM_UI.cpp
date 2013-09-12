#include "AniEditDM_UI.h"
using namespace LSW_ANI_EDIT_UI;

AniEditDM_UI::AniEditDM_UI(QMainWindow *main_win, pAniEditDM data_model):
  data_model(data_model){

  this->file_dialog = pFileDialog(new FileDialog(main_win));
}
	
// interpolate
void AniEditDM_UI::interpolate(){

  file_dialog->warning(data_model != NULL && data_model->interpolate());
}

// update reference sequence
void AniEditDM_UI::updateRefSeq(){
  
  file_dialog->warning( data_model != NULL && data_model->updateRefSeq() );
}

// manager constraint nodes
void AniEditDM_UI::removeAllPosCon(){
  
  if (data_model)
	data_model->removeAllPosCon();
}

// record drag
void AniEditDM_UI::toggleRecord(){
  
  if (data_model){
	const bool record = data_model->isRecordDrag() ? false : true;
	data_model->setRecord(record);
  }
}

// file io
void AniEditDM_UI::saveParitalCon()const{
  
  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveParitalCon(filename));
  }
}
	
void AniEditDM_UI::loadParitalCon(){

  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(data_model != NULL && data_model->loadParitalCon(filename));
	emit update();
  }
}

void AniEditDM_UI::loadConPath(){

  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(data_model != NULL && data_model->loadConPath(filename));
	emit update();
  }
}
	
void AniEditDM_UI::saveOutputMeshes()const{

  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model != NULL && data_model->saveOutputMeshes(filename));
  }
}

void AniEditDM_UI::saveInputMeshes()const{

  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model != NULL && data_model->saveInputMeshes(filename));
  }
}

void AniEditDM_UI::saveCurrentOutputMesh()const{
  
  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveCurrentOutputMesh(filename));
  }
}

void AniEditDM_UI::saveCurrentInputMesh()const{
  
  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveCurrentInputMesh(filename));
  }
}

void AniEditDM_UI::saveCurrentOutputVolMesh(){

  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveCurrentOutputVolMesh(filename));
  }
}

void AniEditDM_UI::saveVolFullU()const{
  
  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveVolFullU(filename));
  }
}

void AniEditDM_UI::saveDragRecord()const{
  
  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveDragRecord(filename));
  }
}

void AniEditDM_UI::loadDragRecord(){

  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->loadDragRecord(filename));
  }
}

void AniEditDM_UI::playRecodDrag(){
  
  file_dialog->warning(data_model != NULL && data_model->playRecodDrag());
}

void AniEditDM_UI::loadCurrentReducedEdits(){
  
  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->loadCurrentReducedEdits(filename));
  }
}

void AniEditDM_UI::saveCurrentReducedEdits()const{
  
  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveCurrentReducedEdits(filename));
  }
}

void AniEditDM_UI::loadAllReducedEdits()const{
  
  const string filename = file_dialog->load();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->loadAllReducedEdits(filename));
  }
}

void AniEditDM_UI::saveAllReducedEdits()const{
  
  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveAllReducedEdits(filename));
  }
}

void AniEditDM_UI::saveRotModalMatrix()const{
  
  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&
						 data_model->getInterpolator()!= NULL&&
						 data_model->getInterpolator()->save("rotational_modal_matrix",filename));
  }
}

void AniEditDM_UI::printEigenValues()const{
  
  if (data_model && data_model->getInterpolator()){
	data_model->getInterpolator()->print("eigenvalues");
  }
}

void AniEditDM_UI::saveCurrentVolFullU()const{

  const string filename = file_dialog->save();
  if(filename.size() >0){
	file_dialog->warning(data_model!=NULL&&data_model->saveCurrentVolFullU(filename));
  }
}
