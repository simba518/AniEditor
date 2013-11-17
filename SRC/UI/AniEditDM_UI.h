#ifndef _ANIEDITDM_UI_H_
#define _ANIEDITDM_UI_H_

#include <QObject>
#include <QtGui/QMainWindow>
#include <FileDialog.h>
#include "AniEditDM.h"
using namespace QGLVEXT;

namespace ANI_EDIT_UI{

#define ANIEDITDM_SAVE(functionName)									\
  const string fname = file_dialog->save();								\
  if(fname.size() >0){													\
	file_dialog->warning(data_model!=NULL&&data_model->functionName(fname)); \
  }

#define ANIEDITDM_LOAD(functionName)									\
  const string fname = file_dialog->load();								\
  if(fname.size() >0){													\
	file_dialog->warning(data_model!=NULL&&data_model->functionName(fname)); \
	emit update();														\
  }

  class AniEditDM_UI:public QObject{

	Q_OBJECT

  public:
	AniEditDM_UI(QMainWindow *main_win, pAniEditDM data_model):
	  data_model(data_model){
	  this->file_dialog = pFileDialog(new FileDialog(main_win));
	}
	
  public slots: 
	void toggleWarp(){
	  if(data_model){
		data_model->useWarp( !(data_model->isUseWarp()) );
		emit update();
	  }
	}
	void interpolate(){
	  file_dialog->warning(data_model != NULL && data_model->interpolate());
	}
	void removeAllPosCon(){
	  if (data_model) data_model->removeAllPosCon();
	}
	void computeConNodeTrajectory(){
	  if(data_model){
		data_model->computeConNodeTrajectory();
		emit update();
	  }
	}
	void toggleRecord(){
	  if (data_model){
		const bool record = data_model->isRecordDrag() ? false : true;
		data_model->setRecord(record);
	  }
	}
	void playRecodDrag(){
	  file_dialog->warning(data_model != NULL && data_model->playRecodDrag());
	}

	void saveParitalCon()const{ANIEDITDM_SAVE(saveParitalCon);}
	void saveOutputVolMeshesVTK()const{ANIEDITDM_SAVE(saveOutputVolMeshesVTK);}
	void saveOutputMeshes()const{ANIEDITDM_SAVE(saveOutputMeshes);}
	void saveInputMeshes()const{ANIEDITDM_SAVE(saveInputMeshes);}
	void saveCurrentOutputMesh()const{ANIEDITDM_SAVE(saveCurrentOutputMesh);}
	void saveCurrentInputMesh()const{ANIEDITDM_SAVE(saveCurrentInputMesh);}
	void saveCurrentOutputVolMesh(){ANIEDITDM_SAVE(saveCurrentOutputVolMesh);}
	void saveVolFullU()const{ANIEDITDM_SAVE(saveVolFullU);}
	void saveCurrentReducedEdits()const{ANIEDITDM_SAVE(saveCurrentReducedEdits);}
	void saveAllReducedEdits()const{ANIEDITDM_SAVE(saveAllReducedEdits);}
	void saveCurrentVolFullU()const{ANIEDITDM_SAVE(saveCurrentVolFullU);}
	void saveDragRecord()const{ANIEDITDM_SAVE(saveDragRecord);}

	void loadParitalCon(){ANIEDITDM_LOAD(loadParitalCon);}
	void loadConPath(){ANIEDITDM_LOAD(loadConPath);}
	void loadCurrentReducedEdits(){ANIEDITDM_LOAD(loadCurrentReducedEdits);}
	void loadAllReducedEdits(){ANIEDITDM_LOAD(loadAllReducedEdits);}
	void loadDragRecord(){ANIEDITDM_LOAD(loadDragRecord);}

	void printEigenValues()const{
	  if (data_model && data_model->getInterpolator()){
		data_model->getInterpolator()->print("eigenvalues");
	  }
	}

  signals:
	void update();

  private:
	pAniEditDM data_model;
	pFileDialog file_dialog;
  };
  typedef boost::shared_ptr<AniEditDM_UI> pAniEditDM_UI;
  
}//end of namespace

#endif /*_ANIEDITDM_UI_H_*/
