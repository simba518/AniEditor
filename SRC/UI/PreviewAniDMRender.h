#ifndef _PREVIEWANIDMRENDER_H_
#define _PREVIEWANIDMRENDER_H_

#include <boost/shared_ptr.hpp>
#include <SelfRenderEle.h>
#include <TextForRender.h>
#include <QObject>
#include <AniCtrl.h>
#include <QGLViewerExt.h>
#include <MeshRender.h>
#include "PreviewAniDM.h"
using namespace QGLVEXT;

namespace ANI_EDIT_UI{
  
  /// @class PreviewAniDMRender render the data of the class PreviewAniDM.
  class PreviewAniDMRender: public SelfRenderEle, public TextForRender{
	
  public:
	PreviewAniDMRender(pPreviewAniDM data_model):data_model(data_model){

	  input_obj_fAmbient[0] = 0.2f;
	  input_obj_fAmbient[1] = 0.2f;
	  input_obj_fAmbient[2] = 0.2f;
	  input_obj_fAmbient[3] = 1.0f;

	  input_obj_fDiffuse[0] = 0.8f;
	  input_obj_fDiffuse[1] = 0.8f;
	  input_obj_fDiffuse[2] = 0.8f;
	  input_obj_fDiffuse[3] = 1.0f;
	}
	void render(){
	  updateFrameNumText();
	  this->draw();
	}
	void draw()const{
	  drawInputObj();
	  drawOutputObj();
	}

  protected:
	void updateFrameNumText(){
  
	  if (data_model != NULL){
	
		const int x = 60;
		const int y = 70;
		const string pre = "Preview";
		this->update(pre,x,y);

		const string frame_num = boost::lexical_cast<string>(data_model->currentFrameNum());
		const string T = boost::lexical_cast<string>( data_model->totalFrameNum() -1 );
		const string text = string("Timestep  ") + frame_num + string(" / ")+T;
		this->update(text,x,y+45);
	  }
	}
	void drawInputObj()const{
	  if (data_model != NULL){
		pObjmesh_const obj_mesh= data_model->getInputObjMesh();
		if (obj_mesh != NULL){
		  UTILITY::draw(*obj_mesh);
		}
	  }
	}
	void drawOutputObj()const{

	  if (data_model != NULL){
		pObjmesh_const obj_mesh= data_model->getOutputObjMesh();
		if (obj_mesh != NULL){
		  UTILITY::draw(*obj_mesh);
		}
	  }
	}
	
  private:
	pPreviewAniDM data_model;
	float input_obj_fAmbient[4];
	float input_obj_fDiffuse[4];
  };
  typedef boost::shared_ptr<PreviewAniDMRender> pPreviewAniDMRender;

  /// @class PreviewAniDMRenderCtrl control the rendering of PreviewAniDM.
  class PreviewAniDMRenderCtrl: public QObject{

	Q_OBJECT

  public:
	PreviewAniDMRenderCtrl(pQGLViewerExt viewer,pAniEditDM dm):viewer(viewer){

	  pPreviewAniDM preview_data_model = pPreviewAniDM(new PreviewAniDM(dm));
	  render = pPreviewAniDMRender(new PreviewAniDMRender(preview_data_model) );
	  if(viewer != NULL){

	  	viewer->addSelfRenderEle(render);
	  	viewer->addTextForRender(render);
	  	preview_ani_ctrl = pAniCtrl(new AniCtrl(preview_data_model,viewer));
	  	preview_ani_ctrl->setPlayCircle(true);
	  }
	}
	
  private:
	pPreviewAniDMRender render;
	pQGLViewerExt viewer;
	pAniCtrl preview_ani_ctrl;
  };
  typedef boost::shared_ptr<PreviewAniDMRenderCtrl> pPreviewAniDMRenderCtrl;
  
}//end of namespace

#endif /*_PREVIEWANIDMRENDER_H_*/
