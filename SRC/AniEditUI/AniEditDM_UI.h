#ifndef _ANIEDITDM_UI_H_
#define _ANIEDITDM_UI_H_

#include <QObject>
#include <QtGui/QMainWindow>
#include <FileDialog.h>
#include "AniEditDM.h"
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{
  
  /**
   * @class AniEditDM_UI provide singnals and slots for AniEditDM.
   * 
   */
  class AniEditDM_UI:public QObject{

	Q_OBJECT

  public:
	AniEditDM_UI(QMainWindow *main_win, pAniEditDM data_model);
	
  public slots: 
	// warping
	void toggleWarp(){

	  if(data_model){
		data_model->useWarp( !(data_model->isUseWarp()) );
		emit update();
	  }
	}

	// interpolate
	void interpolate();

	// update reference sequence
	void updateRefSeq();

	// manager constraint nodes
	void removeAllPosCon();
	void computeConNodeTrajectory(){
	  if(data_model){
		data_model->computeConNodeTrajectory();
		emit update();
	  }
	}

	// record drag
	void toggleRecord();
	void playRecodDrag();

	// file io
	void saveParitalCon()const;
	void loadParitalCon();
	void loadConPath();
	void saveOutputMeshes()const;
	void saveInputMeshes()const;
	void saveCurrentOutputMesh()const;
	void saveCurrentInputMesh()const;
	void saveCurrentOutputVolMesh();
	void saveVolFullU()const;
	void loadCurrentReducedEdits();
	void saveCurrentReducedEdits()const;
	void loadAllReducedEdits()const;
	void saveAllReducedEdits()const;	
	void saveRotModalMatrix()const;
	void saveCurrentVolFullU()const;

	void saveDragRecord()const;
	void loadDragRecord();

	void printEigenValues()const;

  signals:
	void update();

  private:
	pAniEditDM data_model;
	pFileDialog file_dialog;
  };
  
  typedef boost::shared_ptr<AniEditDM_UI> pAniEditDM_UI;
  
}//end of namespace

#endif /*_ANIEDITDM_UI_H_*/
