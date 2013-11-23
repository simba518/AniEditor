#ifndef _ANIEDITMAINWIN_H_
#define _ANIEDITMAINWIN_H_

#include <QtGui/QMainWindow>
#include <FileDialog.h>
#include <QGLViewerExt.h>
#include <AniCtrl.h>
#include <AniSliderBar.h>
using namespace QGLVEXT;

#include <SRC/UI/ui_aniedit.h>
#include <SRC/UI/ui_about.h>
#include "VolObjCtrl.h"
#include "AniEditDM_UI.h"
#include "AniEditDMSelection.h"
#include "AniEditDMRender.h"
#include "AniEditDMDragCtrl.h"
#include "PreviewAniDMRender.h"

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

	// animation 
	pAniCtrl ani_ctrl;

	// data model and operations
	pTetMeshEmbeding p_VolObjMesh;
	UTILITY::pVolObjCtrl p_VolObjCtrl;
	pAniEditDM p_AniEditDM;
	pAniDataModel p_animation;
	pAniEditDM_UI p_AniEditDM_UI;
	pAniEditDMRenderCtrl p_AniEditDMRenderCtrl;
	pAniEditDMDragCtrl p_AniEditDMDragCtrl;
	pAniEditDMSelCtrl p_AniEditDMSelCtrl;

	// preview animation
	pPreviewAniDMRenderCtrl p_PreviewAniDMRenderCtrl;
	pQGLViewerExt preview_viewer;

	// init file
	string init_filename; //< recorded for reloading.
	string main_win_initfile; //< the initfile control the appearance of the
							  //window.


  };
  
}//end of namespace

#endif /*_ANIEDITMAINWIN_H_*/
