#ifndef _KEYFRAMEDMRENDER_H_
#define _KEYFRAMEDMRENDER_H_

#include <QObject>
#include <SelfRenderEle.h>
#include <QGLViewerExt.h>
#include <MeshRender.h>
#include <VolNodeGroupRender.h>
#include "KeyframeDM.h"
using namespace QGLVEXT;
using namespace UTILITY;

namespace ANI_EDIT_UI{

  class KeyframeDMRender: public SelfRenderEle{
	
  public:
	KeyframeDMRender(pKeyframeDM dm):data_model(dm){
	  const double con_group_color[4] = {0.6f, 0.0f, 0.0f, 1.0f};
	  node_group_render.setColor(con_group_color);
	  showAllKeyframes = false;
	}
	void draw()const{

	  double x=0,y=0,z=0;
	  if (data_model){
		data_model->getSceneTranslate(x,y,z);
		glTranslated(x,y,z);
	  }
	  if (showAllKeyframes&&data_model){
		const vector<VectorXd> &u = data_model->getAllKeyU();
		for (int i = 0; i < u.size(); ++i){
		  UTILITY::draw(data_model->getObjMesh(i),keyframe_mtl);
		  // drawVolMesh(u[i]);
		}
	  }else if (data_model&&data_model->isKeyframe()>=0) {
		UTILITY::draw(data_model->getObjMesh(),keyframe_mtl);
		// drawVolMesh(data_model->getKeyU());
		// pPartialConstraints_const par_con = data_model->getPartialCon();
		// assert(par_con);
		// const VectorXd &vol_u = data_model->getKeyU();
		// pTetMesh_const tetmesh = data_model->getVolMesh();
		// node_group_render.draw(tetmesh,par_con->getConNodesSet(),&vol_u[0],DRAW_POINT);
	  }
	  glTranslated(-x,-y,-z);
	}
	void toggleShowAllKeyframes(){
	  showAllKeyframes = showAllKeyframes ? false:true;
	}
	
  protected:
	void drawVolMesh(const VectorXd &vol_u)const{

	  assert(data_model);
	  pTetMesh_const tetmesh = data_model->getVolMesh();
	  assert(tetmesh);
	  assert_eq(vol_u.size(),tetmesh->nodes().size()*3);
	  assert_gt(vol_u.size(),0);
	  UTILITY::draw(tetmesh,&(vol_u[0]));
	}

  private:
	UTILITY::VolNodeGroupRender node_group_render;
	bool showAllKeyframes;
	pKeyframeDM data_model;
	ObjMtl keyframe_mtl;
  };
  typedef boost::shared_ptr<KeyframeDMRender> pKeyframeDMRender;

  /// @class AniEditDMRenderCtrl control the rendering of AniEditDM class.
  class KeyframeDMRenderCtrl: public QObject{
	
	Q_OBJECT
	
  public:
	KeyframeDMRenderCtrl(pQGLViewerExt viewer,pKeyframeDM dm):viewer(viewer){
  
	  model_render=pKeyframeDMRender(new KeyframeDMRender(dm));
	  if(viewer != NULL) viewer->addSelfRenderEle(model_render);
	}
	
  public slots:
	void toggleShow(){
	  if (viewer != NULL) viewer->toggleRemoveAddSelfRenderEle(model_render);
	}
	void toggleShowAllKeyframes(){
	  model_render->toggleShowAllKeyframes();
	}

  private:
	pKeyframeDMRender model_render;
	pQGLViewerExt viewer;
  };
  typedef boost::shared_ptr<KeyframeDMRenderCtrl> pKeyframeDMRenderCtrl;

}
#endif /* _KEYFRAMEDMRENDER_H_ */
