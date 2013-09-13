#ifndef _ANIEDITDMDRAGCTRL_H_
#define _ANIEDITDMDRAGCTRL_H_

#include <boost/shared_ptr.hpp>
#include <QObject>
#include <DragCtrl.h>
#include "AniEditDMConNodeDrag.h"
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{
  
  /**
   * @class AniEditDMDragCtrl
   * 
   */
  class AniEditDMDragCtrl:public QObject{

	Q_OBJECT
	
  public: 
	AniEditDMDragCtrl(pQGLViewerExt view,pAniEditDM dm,bool show_circle=false){

	  con_nodes_drag = pAniEditDMConNodeDrag(new AniEditDMConNodeDrag(view,dm));
	  drag_ctrl = pDragCtrl( new DragCtrl(view,con_nodes_drag) );
	  drag_ctrl->setObserver(con_nodes_drag);
	  drag_ctrl->setDragHook(con_nodes_drag);
	  if (view != NULL){
		con_nodes_drag->setRenderPriority(FIRST_RENDER);
		view->addSelfRenderEle(con_nodes_drag);
	  }
	  showMouseCircle(show_circle);
	  showDragCircle(show_circle);
	}
	void showMouseCircle(bool show){
	  con_nodes_drag->showMouseCircle(show);
	}
	void showDragCircle(bool show){
	  con_nodes_drag->showDragCircle(show);
	}

  public slots: 
	// rendering
	void toggleShowMouseCircle(){
	  con_nodes_drag->toggleShowMouseCircle();
	}
	void toggleShowDragCircle(){
	  con_nodes_drag->toggleShowDragCircle();
	}
	
  private:
	pAniEditDMConNodeDrag con_nodes_drag;
	pDragCtrl drag_ctrl;
  };
  
  typedef boost::shared_ptr<AniEditDMDragCtrl> pAniEditDMDragCtrl;
  
}//end of namespace

#endif /*_ANIEDITDMDRAGCTRL_H_*/
