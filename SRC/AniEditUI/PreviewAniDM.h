#ifndef _PREVIEWANIDM_H_
#define _PREVIEWANIDM_H_

#include <boost/shared_ptr.hpp>
#include <AniDataModel.h>
#include "AniEditDM.h"

namespace LSW_ANI_EDIT_UI{
  
  /**
   * @class PreviewAniDM animation data model for preview of the edited
   * animation.
   * 
   */
  class PreviewAniDM: public QGLVEXT::AniDataModel{
	
  public:
	PreviewAniDM(pAniEditDM data_model):data_model(data_model){}
	int totalFrameNum()const{
	  if (data_model != NULL){
		return data_model->totalFrameNum();
	  }
	  return 0;
	}

	pObjRenderMesh_const getInputObjMesh()const{
	  assert (data_model != NULL);
	  return data_model->getInputObjMesh(this->currentFrameNum());
	}
	pObjRenderMesh_const getOutputObjMesh()const{
	  assert (data_model != NULL);
	  return data_model->getOutputObjMesh(this->currentFrameNum());
	}
	
  private:
	pAniEditDM data_model;
  };
  
  typedef boost::shared_ptr<PreviewAniDM> pPreviewAniDM;
  
}//end of namespace

#endif /*_PREVIEWANIDM_H_*/
