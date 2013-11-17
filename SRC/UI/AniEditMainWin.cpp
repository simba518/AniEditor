#include <boost/filesystem.hpp>
#include <Log.h>
#include <JsonFilePaser.h>
#include "AniEditMainWin.h"
using namespace UTILITY;
using namespace ANI_EDIT_UI;

AniEditMainWin::AniEditMainWin(const string w_inif,QWidget *p,Qt::WFlags flags):
  QMainWindow(p,flags),main_win_initfile(w_inif){
  initComponents(p,flags);
  createConnections();
  paserCommandLine();
}

void AniEditMainWin::initComponents(QWidget *parent,Qt::WFlags flags){

  initViewers();

  file_dialog = pFileDialog(new FileDialog(this));

  p_animation = pAniDataModel(new AniDataModel());
  p_VolObjMesh = pTetMeshEmbeding(new TetMeshEmbeding());
  p_VolObjCtrl = pVolObjCtrl(new VolObjCtrl(this, p_VolObjMesh));
  p_AniEditDM = pAniEditDM(new AniEditDM(p_VolObjMesh,p_animation));
  p_AniEditDM_UI = pAniEditDM_UI(new AniEditDM_UI(this, p_AniEditDM));
  p_AniEditDMRenderCtrl = pAniEditDMRenderCtrl(new AniEditDMRenderCtrl(viewer, p_AniEditDM));
  p_AniEditDMDragCtrl = pAniEditDMDragCtrl(new AniEditDMDragCtrl(viewer, p_AniEditDM));
  p_AniEditDMSelCtrl = pAniEditDMSelCtrl(new AniEditDMSelCtrl(viewer, p_AniEditDM));
  ani_ctrl = pAniCtrl(new AniCtrl(p_animation,viewer));
  m_mainwindow.ani_slider->setAniCtrl(ani_ctrl);

  p_PreviewAniDMRenderCtrl = pPreviewAniDMRenderCtrl(new PreviewAniDMRenderCtrl(preview_viewer, p_AniEditDM));
  preview_viewer->hide();
}

void AniEditMainWin::initViewers(){

  m_mainwindow.setupUi(this);
  viewer = m_mainwindow.left_view;
  preview_viewer = m_mainwindow.right_view;
  about_Dialog = new QDialog(this);
  m_about.setupUi(about_Dialog);
  viewer->setAnimationPeriod(1); 
  preview_viewer->setAnimationPeriod(30); 
}

void AniEditMainWin::createConnections(){

  // windows
  connect(m_mainwindow.actionResetWindows, SIGNAL(triggered()), this, SLOT(resetWindows()));
  connect(m_mainwindow.actionPreview, SIGNAL(toggled(bool)), preview_viewer, SLOT(setVisible(bool)));

  // about
  connect(m_mainwindow.actionAbout, SIGNAL(triggered()), about_Dialog, SLOT(show()));

  // load init file
  connect(m_mainwindow.actionLoadInitFile, SIGNAL(triggered()), this, SLOT(loadInitFile()));
  connect(m_mainwindow.actionReloadInitFile, SIGNAL(triggered()), this, SLOT(reloadInitfile()));
  connect(m_mainwindow.actionLoadStateFile, SIGNAL(triggered()), viewer, SLOT(loadStateFile()));
  connect(m_mainwindow.actionSaveStateFile, SIGNAL(triggered()), viewer, SLOT(saveStateFile()));

  // update viewer
  connect(p_VolObjCtrl.get(), SIGNAL(resetSceneMsg(double ,double ,double ,double ,double ,double )),
  		  viewer, SLOT(resetSceneBoundBox(double ,double ,double ,double ,double ,double )));
  connect(p_VolObjCtrl.get(), SIGNAL(resetSceneMsg(double ,double ,double ,double ,double ,double )),
  		  preview_viewer, SLOT(resetSceneBoundBox(double ,double ,double ,double ,double ,double )));
  connect(m_mainwindow.actionResetViewer, SIGNAL(triggered()), p_VolObjCtrl.get(), SLOT(sendMsgResetScene()));
  // connect(m_mainwindow.actionPhoneShading, SIGNAL(triggered()), p_VolObjCtrl.get(), SLOT(togglePhoneShading()));
  
  connect(p_VolObjCtrl.get(), SIGNAL(update()), viewer, SLOT(update()));
  connect(p_VolObjCtrl.get(), SIGNAL(update()), preview_viewer, SLOT(update()));

  // viewer: animation 
  connect(m_mainwindow.actionPlay, SIGNAL(triggered()), viewer, SLOT(startAnimation()));
  connect(m_mainwindow.actionStop, SIGNAL(triggered()), viewer, SLOT(stopAnimation()));
  connect(m_mainwindow.actionPause, SIGNAL(triggered()), viewer, SLOT(pauseAnimation()));
  connect(m_mainwindow.actionStepByStep, SIGNAL(triggered()), viewer, SLOT(stepByStep()));
  connect(m_mainwindow.actionBackward, SIGNAL(triggered()), viewer, SLOT(backward()));
  connect(m_mainwindow.actionForward, SIGNAL(triggered()), viewer, SLOT(forward()));
  connect(m_mainwindow.actionFastforward, SIGNAL(triggered()), viewer, SLOT(fastForward()));
  connect(m_mainwindow.actionPushViewerStatus, SIGNAL(triggered()), viewer, SLOT(pushStatus()));
  connect(m_mainwindow.actionRestoreViewerStatus, SIGNAL(triggered()), viewer, SLOT(restoreStatus()));
  connect(m_mainwindow.actionPopViewerStatus, SIGNAL(triggered()), viewer, SLOT(popStatus()));

  // QtBaseAniControl
  connect(m_mainwindow.actionPlayCircle, SIGNAL(triggered()), ani_ctrl.get(), SLOT(togglePlayCircle()));

  // preview viewer: animation 
  connect(m_mainwindow.actionPlay_preview, SIGNAL(triggered()), preview_viewer, SLOT(startAnimation()));
  connect(m_mainwindow.actionPause_preview, SIGNAL(triggered()), preview_viewer, SLOT(pauseAnimation()));

  // AniEditDM_UI
  connect(m_mainwindow.actionInterpolate, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(interpolate()));
  connect(m_mainwindow.actionRemoveAllConNodes, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(removeAllPosCon()));
  connect(m_mainwindow.actionSavePartialCon, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveParitalCon()));
  connect(m_mainwindow.actionLoadPartialCon, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(loadParitalCon()));
  connect(m_mainwindow.actionLoadConPath, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(loadConPath()));
  connect(m_mainwindow.actionSaveOutputVolMeshesVTK, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveOutputVolMeshesVTK()));
  connect(m_mainwindow.actionSaveOutputMeshes, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveOutputMeshes()));
  connect(m_mainwindow.actionSaveInputMeshes, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveInputMeshes()));
  connect(m_mainwindow.actionSaveCurrentOutputMesh, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveCurrentOutputMesh()));
  connect(m_mainwindow.actionSaveCurrentInputMesh, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveCurrentInputMesh()));
  connect(m_mainwindow.actionSaveCurrentOutputVolMesh, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveCurrentOutputVolMesh()));
  connect(m_mainwindow.actionSaveVolFullU, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveVolFullU()));
  connect(m_mainwindow.actionUseWarp, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(toggleWarp()));
  connect(p_AniEditDM_UI.get(), SIGNAL(update()), viewer, SLOT(update()));

  connect(m_mainwindow.actionSaveCuReducedEdit, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveCurrentReducedEdits()));
  connect(m_mainwindow.actionSaveAllReducedEdit, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveAllReducedEdits()));
  connect(m_mainwindow.actionLoadAllReducedEdit, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(loadAllReducedEdits()));
  connect(m_mainwindow.actionSaveCurrentVolFullU, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveCurrentVolFullU()));
  connect(m_mainwindow.actionLoadCuReducedEdit, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(loadCurrentReducedEdits()));
  // connect(m_mainwindow.actionSaveRotModalMat, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveRotModalMatrix()));
  connect(m_mainwindow.actionPrintEigenvalues, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(printEigenValues()));
  connect(m_mainwindow.actionUpdateConTrajectory, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(computeConNodeTrajectory()));
  

  // AniEditDM_UI: record drag operation
  connect(m_mainwindow.actionRecordDragOp, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(toggleRecord()));
  connect(m_mainwindow.actionPlayDragOp, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(playRecodDrag()));
  connect(m_mainwindow.actionSaveDragOp, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(saveDragRecord()));
  connect(m_mainwindow.actionLoadDragOp, SIGNAL(triggered()), p_AniEditDM_UI.get(), SLOT(loadDragRecord()));

  // AniEditDMRenderCtrl
  connect(m_mainwindow.actionShowInputObj, SIGNAL(triggered()), p_AniEditDMRenderCtrl.get(), SLOT(toggleShowInputObj()));
  connect(m_mainwindow.actionShowInputVol, SIGNAL(triggered()), p_AniEditDMRenderCtrl.get(), SLOT(toggleShowInputVol()));
  connect(m_mainwindow.actionShowOutputObj, SIGNAL(triggered()), p_AniEditDMRenderCtrl.get(), SLOT(toggleShowOutputObj()));
  connect(m_mainwindow.actionShowOutputVol, SIGNAL(triggered()), p_AniEditDMRenderCtrl.get(), SLOT(toggleShowOutputVol()));
  //connect(m_mainwindow.actionShowKeyframe, SIGNAL(triggered()), p_AniEditDMRenderCtrl.get(), SLOT(toggleShowKeyframeObj()));
  connect(m_mainwindow.actionShowConNodes, SIGNAL(triggered()), p_AniEditDMRenderCtrl.get(), SLOT(toggleShowConNodes()));
  connect(m_mainwindow.actionShowConPath, SIGNAL(triggered()), p_AniEditDMRenderCtrl.get(), SLOT(toggleShowConPath()));

  // AniEditDMDragCtrl
  // connect(m_mainwindow.actionEnableDrag, SIGNAL(triggered()), p_AniEditDMDragCtrl.get(), SLOT(toggleDrag()));
  connect(m_mainwindow.actionMouseCircle, SIGNAL(triggered()), p_AniEditDMDragCtrl.get(), SLOT(toggleShowMouseCircle()));
  connect(m_mainwindow.actionConCircle, SIGNAL(triggered()), p_AniEditDMDragCtrl.get(), SLOT(toggleShowDragCircle()));

  // AniEditDMSelCtrl 
  connect(m_mainwindow.actionPrintSelNodes, SIGNAL(triggered()), p_AniEditDMSelCtrl.get(), SLOT(togglePrintSelEle()));
  connect(m_mainwindow.actionEnableSelect, SIGNAL(triggered()), p_AniEditDMSelCtrl.get(), SLOT(toggleSelOp()));


  // draw lights
  connect(m_mainwindow.actionLights, SIGNAL(triggered()), viewer, SLOT(toggleDrawLights()));  
  connect(m_mainwindow.actionLights, SIGNAL(triggered()), preview_viewer, SLOT(toggleDrawLights()));  
}

void AniEditMainWin::paserCommandLine(){

  QStringList cmdline_args = QCoreApplication::arguments();
  if(cmdline_args.size() >= 2 ){

  	string init_filename = cmdline_args.at(1).toStdString();
	bool succ = false;
  	if(boost::filesystem::exists(init_filename)){
  	  succ = loadInitFile(init_filename);
  	  if(succ){
		INFO_LOG("success to load initfile from  " << init_filename);
  	  }else{
		ERROR_LOG("failed to load initfile from  "<< init_filename);
  	  }	  
  	}else{
	  ERROR_LOG("file " << init_filename <<" is not exist!" );
  	}
	autoRun(init_filename);
  }
}

bool AniEditMainWin::loadInitFile(const string init_filename){

  TRACE_FUN();
  
  bool succ = false;
  if (p_VolObjCtrl != NULL){
	succ = p_VolObjCtrl->initialize(init_filename);
  }
  if (succ && p_AniEditDM != NULL){
	succ = p_AniEditDM->initialize(init_filename);
  }
  if (succ){
	succ = paserInitFile(init_filename);
  }
  if (succ){
	p_animation->reset();
  }
  p_AniEditDM->print();
  return succ;
}

void AniEditMainWin::loadInitFile(){
  
  const string filename = file_dialog->load("ini");
  if(filename.size() >0){
	assert (file_dialog != NULL);
	file_dialog->warning(loadInitFile(filename));
  }
}

void AniEditMainWin::reloadInitfile(){

  if(init_filename.size() > 0){
	file_dialog->warning(loadInitFile(init_filename),"failed to reload the initfile!");
  }
}

bool AniEditMainWin::paserInitFile(const string init_filename){
  
  return true;
}

void AniEditMainWin::resetWindows(){

  const int default_h = 800;
  const int default_w = 720;
  if (viewer) 
	viewer->resize(default_w,default_h);
  if (preview_viewer) 
	preview_viewer->resize(default_w,default_h);
  this->resize(default_w*2,default_h);
}

void AniEditMainWin::autoRun(const string init_filename){
  
  UTILITY::JsonFilePaser inf;
  if(inf.open(init_filename)){
	
	bool auto_run = false;
	bool auto_exit = false;
	string interp_seq_save_file;

	if(inf.read("auto_run",auto_run) && auto_run){
	  
	  if (p_AniEditDM){
		const bool succ = p_AniEditDM->playRecodDrag();
		ERROR_LOG_COND("failed to perform interpolation ",succ);
	  }
	  if(inf.read("auto_exit",auto_exit) && auto_exit){
		INFO_LOG("auto exit, controlled by the json node 'auto_exit' of the file " << init_filename );
		exit(0);
	  }
	}
  }

}
