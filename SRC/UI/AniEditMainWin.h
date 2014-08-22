#ifndef _ANIEDITMAINWIN_H_
#define _ANIEDITMAINWIN_H_

#include <QtGui/QMainWindow>
#include <FileDialog.h>
#include <QGLViewerExt.h>
#include <ConTrackBall.h>
#include <AniCtrl.h>
#include <AniSliderBar.h>
using namespace QGLVEXT;

#ifdef WIN32
#include <ui_aniedit.h>
#include <ui_about.h>
#else
#include <SRC/UI/ui_aniedit.h>
#include <SRC/UI/ui_about.h>
#endif

#include "VolObjCtrl.h"
#include "AniEditDM_UI.h"
#include "AniEditDMSelection.h"
#include "AniEditDMRender.h"
#include "AniEditDMConNodeDrag.h"
#include "PreviewAniDMRender.h"
#include "RecordAndReplayWrapper.h"
#include "KeyframeDM.h"
#include "KeyframeDMRender.h"
#include "KeyframeDMSelection.h"
#include "ManipulateOP.h"
using namespace UTILITY;

namespace ANI_EDIT_UI{
  
  /**
   * @class AniEditMainWin main window for animation edit.
   * 
   */
  class AniEditMainWin: public QMainWindow{

	Q_OBJECT
	
  public:
	AniEditMainWin(const string win_inif,QWidget *parent=0,Qt::WFlags flags=0);
	bool loadInitFile(const string init_filename);
	
  protected:
	void initComponents(QWidget *parent,Qt::WFlags flags);
	void createConnections();
	void paserCommandLine();
	bool paserInitFile(const string init_filename);
	void initViewers();

  protected slots:
	void loadInitFile();
	void reloadInitfile();
	void resetWindows();

  private:
	// view
	Ui::MainWindow m_mainwindow;
	Ui::About m_about;
	pQGLViewerExt viewer;
	pFileDialog file_dialog;
	QDialog *about_Dialog;
	pConTrackBall p_ConTrackBall;

	// animation 
	pAniCtrl ani_ctrl;

	// data model and operations
	pTetMeshEmbeding p_VolObjMesh;
	UTILITY::pVolObjCtrl p_VolObjCtrl;
	pAniEditDM p_AniEditDM;
	pLocalframeManipulatoion manipulation;
	pLocalframeManipulatoionCtrl manipulation_ctrl;
	pAniDataModel p_animation;
	pAniEditDM_UI p_AniEditDM_UI;
	pAniEditDMRenderCtrl p_AniEditDMRenderCtrl;
	pAniEditDMDragCtrl p_AniEditDMDragCtrl;
	pAniEditDMSelCtrl p_AniEditDMSelCtrl;
	pKeyframeDMRenderCtrl p_KeyframeDMRenderCtrl;
	pKeyframeDM_UI p_KeyframeDM_UI;
	pKeyframeDMSelectionCtrl p_KeyframeDMSelectionCtrl;

	// preview animation
	pPreviewAniDMRenderCtrl p_PreviewAniDMRenderCtrl;
	pQGLViewerExt preview_viewer;

	//record and replay
	UTILITY::pRecordAndReplayWrapper p_RecordAndReplayWrapper;

	// init file
	string init_filename; //< recorded for reloading.
	string main_win_initfile; //< the initfile control the appearance of the
							  //window.


  };
  
}//end of namespace

#endif /*_ANIEDITMAINWIN_H_*/
