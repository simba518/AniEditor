#ifndef _PREVIEWANIDMRENDER_H_
#define _PREVIEWANIDMRENDER_H_

#include <boost/shared_ptr.hpp>
#include <SelfRenderEle.h>
#include "PreviewAniDM.h"
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{
  
  /**
   * @class PreviewAniDMRender render the data of the class PreviewAniDM.
   * 
   */
  class PreviewAniDMRender: public SelfRenderEle{
	
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
	
		/// @todo
		// const int x = 60;
		// const int y = 70;
		// const string pre = "Preview";
		// this->updateText(pre,x,y);

		// const string frame_num = boost::lexical_cast<string>(data_model->currentFrameNum());
		// const string T = boost::lexical_cast<string>( data_model->totalFrameNum() -1 );
		// const string text = string("Timestep  ") + frame_num + string(" / ")+T;
		// this->updateText(text,x,y+45);
	  }
	}
	void drawInputObj()const{
	  if (data_model != NULL){
		pObjRenderMesh_const obj_mesh= data_model->getInputObjMesh();
		if (obj_mesh != NULL){
		  obj_mesh->draw(0.5,0.5,0.5);
		}
	  }
	}
	void drawOutputObj()const{

	  if (data_model != NULL){
		pObjRenderMesh_const obj_mesh= data_model->getOutputObjMesh();
		if (obj_mesh != NULL){
		  obj_mesh->draw();
		}
	  }
	}
	
  private:
	pPreviewAniDM data_model;
	float input_obj_fAmbient[4];
	float input_obj_fDiffuse[4];
  };
  
  typedef boost::shared_ptr<PreviewAniDMRender> pPreviewAniDMRender;
  
}//end of namespace

#endif /*_PREVIEWANIDMRENDER_H_*/
