#ifndef _ANIEDITDMCONNODEDRAG_H_
#define _ANIEDITDMCONNODEDRAG_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <QGLViewerExt.h>
#include <SelectCtrl.h>
#include <DragCtrl.h>
#include <Selectable.h>
#include "VolNodeGroupRender.h"
#include "AniEditDM.h"
using namespace Eigen;
using namespace UTILITY;
using namespace QGLVEXT;

namespace ANI_EDIT_UI{

  /**
   * @class AniEditDMConNodeDrag provide interfaces for draging the constrained
   * nodes for AniEditDM data model.
   * 
   */
  class AniEditDMConNodeDrag:public Selectable, public DragObserver, 
							 public SelfRenderEle, public DragHook{

  public:
	AniEditDMConNodeDrag(pQGLViewerExt viewer, pAniEditDM data_model);
	
	// drag
	int totalEleNum ()const{
	  return data_model != NULL ? data_model->getConNodes().size() : 0;
	}
	void drawWithNames ()const;
	void prepareSelection(){
	  draggedGroupId = -1;
	}
	void getDragedPoint(double point[3])const;
	
	// observer
	void selectDragEle(int sel_group_id);
	void startDrag (double x,double y,double z);
	void dragTo (double x,double y,double z);
	void stopDrag (double x,double y,double z){
	  draggedGroupId = -1;
	}

	// render
	void render();
	void draw()const;

  protected:
	bool hasDraggedGroup()const{
	  return draggedGroupId >= 0;
	}

  private:
	pQGLViewerExt viewer;
	pAniEditDM data_model;
	int draggedGroupId;

	VolNodeGroupRender con_node_render;
	Vector3d initial_dragged_point; // start position of the dragged point.
	Matrix<double,3,-1> initial_displacement; // initial displacements of the dragged group.
  };
  typedef boost::shared_ptr<AniEditDMConNodeDrag> pAniEditDMConNodeDrag;

  class AniEditDMDragCtrl:public QObject{

	Q_OBJECT
	
  public: 
	AniEditDMDragCtrl(pQGLViewerExt view,pAniEditDM dm){
	  con_nodes_drag = pAniEditDMConNodeDrag(new AniEditDMConNodeDrag(view,dm));
	  drag_ctrl = pDragCtrl( new DragCtrl(view,con_nodes_drag) );
	  drag_ctrl->setObserver(con_nodes_drag);
	  drag_ctrl->setDragHook(con_nodes_drag);
	  drag_ctrl->setKeyMouse(Qt::ControlModifier,Qt::LeftButton);
	  if (view != NULL){
		con_nodes_drag->setRenderPriority(FIRST_RENDER);
		view->addSelfRenderEle(con_nodes_drag);
	  }
	}

  private:
	pAniEditDMConNodeDrag con_nodes_drag;
	pDragCtrl drag_ctrl;
  };
  typedef boost::shared_ptr<AniEditDMDragCtrl> pAniEditDMDragCtrl;
  
}//end of namespace

#endif /*_ANIEDITDMCONNODEDRAG_H_*/
