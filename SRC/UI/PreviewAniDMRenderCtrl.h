#ifndef _PREVIEWANIDMRENDERCTRL_H_
#define _PREVIEWANIDMRENDERCTRL_H_

#include <boost/shared_ptr.hpp>
#include <QObject>
#include <AniCtrl.h>
#include <QGLViewerExt.h>
#include "PreviewAniDMRender.h"
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{
  
  /**
   * @class PreviewAniDMRenderCtrl control the rendering of PreviewAniDM.
   * 
   */
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

#endif /*_PREVIEWANIDMRENDERCTRL_H_*/
