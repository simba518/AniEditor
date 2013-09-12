#ifndef _ANIEDITDMRENDERCTRL_H_
#define _ANIEDITDMRENDERCTRL_H_

#include <boost/shared_ptr.hpp>
#include <QObject>
#include <QGLViewerExt.h>
#include "AniEditDMRender.h"
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{
  
  /**
   * @class AniEditDMRenderCtrl control the rendering of AniEditDM class.
   * 
   */
  class AniEditDMRenderCtrl: public QObject{
	
	Q_OBJECT
	
  public:
	AniEditDMRenderCtrl(pQGLViewerExt viewer,pAniEditDM dm, const string title = "Editing"):viewer(viewer){
  
	  model_render=pAniEditDMRender(new AniEditDMRender(dm,OUTPUT_OBJ|CON_NODES|CON_TRAJECTORY,title));
	  if(viewer != NULL){
		viewer->addSelfRenderEle(model_render);
		viewer->addTextForRender(model_render);
	  }
	}
	
  public slots:
	void toggleShowInputObj(){
	  if (model_render != NULL)
		model_render->toggleRenderType(INPUT_OBJ);
	  if (viewer != NULL) viewer->update();
	}
	void toggleShowInputVol(){
	  if (model_render != NULL)
		model_render->toggleRenderType(INPUT_VOL);
	  if (viewer != NULL) viewer->update();
	}
	void toggleShowOutputObj(){
	  if (model_render != NULL)
		model_render->toggleRenderType(OUTPUT_OBJ);
	  if (viewer != NULL) viewer->update();
	}
	void toggleShowOutputVol(){
	  if (model_render != NULL)
		model_render->toggleRenderType(OUTPUT_VOL);
	  if (viewer != NULL) viewer->update();
	}
	void toggleShowConNodes(){
	  if (model_render != NULL)
		model_render->toggleRenderType(CON_NODES);
	  if (viewer != NULL) viewer->update();
	}
	void toggleShowConPath(){
	  if (model_render != NULL)
		model_render->toggleRenderType(CON_PATH);
	  if (viewer != NULL) viewer->update();
	}

	void showInputSeq(bool show){ showType(show, INPUT_OBJ); }
	void showConNodes(bool show){ showType(show, CON_NODES); }
	void showConPath(bool show){ showType(show, CON_PATH); }

  protected:
	void showType(bool show, ANIEDITDM_RENDER_TYPE type){
	  if(show && model_render != NULL){
		model_render->addRenderType(type);
	  }else if (model_render != NULL){
		model_render->removeRenderType(type);
	  }
	  if (viewer != NULL) viewer->update();
	}
	
  private:
	pAniEditDMRender model_render;
	pQGLViewerExt viewer;
	
  };
  
  typedef boost::shared_ptr<AniEditDMRenderCtrl> pAniEditDMRenderCtrl;
  
}//end of namespace

#endif /*_ANIEDITDMRENDERCTRL_H_*/
