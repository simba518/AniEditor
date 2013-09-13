#ifndef _ANIEDITDMCONNODEDRAG_H_
#define _ANIEDITDMCONNODEDRAG_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <QGLViewerExt.h>
#include <SelectCtrl.h>
#include <DragCtrl.h>
#include <Selectable.h>
#include <BarycenterDraggerVec.h>
#include "VolNodeGroupRender.h"
#include "AniEditDM.h"
using namespace Eigen;
using namespace LSW_SIM;
using namespace LSW_BASE_UI;
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{

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
	void prepareSelection(){
	  if (data_model != NULL&& data_model->currentFrameNum() >=0 ){
		const vector<set<int> > con_nodes = data_model->getConNodes();
		dragger.setConGroups(con_nodes,data_model->getUforConstraint());
	  }
	}
	void drawWithNames ()const;
	void getDragedPoint(double point[3])const;
	
	// observer
	void selectDragEle(int sel_group_id);
	void startDrag(int screen_x, int screen_y){}
	void startDrag (double x,double y,double z){
	  dragger.startDrag(x,y,z);
	}
	void dragTo (double x,double y,double z){
	  dragger.dragTo(x,y,z);
	  updateConPos();
	}
	void stopDrag (double x,double y,double z){
	  dragger.stopDrag();
	  updateConPos();
	}

	// render
	void render();
	void draw()const;
	void showMouseCircle(bool show){
	  draw_mouse_circle = show;
	}
	void showDragCircle(bool show){
	  draw_drag_circle = show;
	}
	void toggleShowMouseCircle(){
	  draw_mouse_circle = draw_mouse_circle ? false:true;
	}
	void toggleShowDragCircle(){
	  draw_drag_circle = draw_drag_circle ? false:true;
	}

  protected:
	void updateConPos();
	const Vector3d getBaryCenter(const set<int>&one_group,const VectorXd&u)const;
	void draw2DCircle(const double screen_width, const double screen_height,
					  const double screen_x, const double screen_y, 
					  const double color[3])const;
	void drawMouseCircle()const;
	void drawDragPointCircle()const;
	
  private:
	pQGLViewerExt viewer;
	pAniEditDM data_model;
	BarycenterDraggerVec dragger;
	VolNodeGroupRender con_node_render;

	Vector3d drag_point; // return the position of the dragged point, which
						  // should be initialized when selectDragEle(int) is
						  // called.
	Vector3d mouse_point; // center of the mouse's circle
	Vector3d draged_point;// center of the draged point's circle
	bool draw_mouse_circle;
	bool draw_drag_circle;
  };
  
  typedef boost::shared_ptr<AniEditDMConNodeDrag> pAniEditDMConNodeDrag;
  
}//end of namespace

#endif /*_ANIEDITDMCONNODEDRAG_H_*/
