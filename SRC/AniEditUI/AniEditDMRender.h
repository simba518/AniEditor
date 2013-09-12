#ifndef _ANIEDITDMRENDER_H_
#define _ANIEDITDMRENDER_H_

#include <boost/shared_ptr.hpp>
#include <SelfRenderEle.h>
#include "VolNodeGroupRender.h" ///@todo
#include "AniEditDM.h"
using namespace QGLVEXT;

namespace LSW_ANI_EDIT_UI{

  enum ANIEDITDM_RENDER_TYPE {NOTHING=0,
							  INPUT_OBJ=1,
							  INPUT_VOL=2,
							  OUTPUT_OBJ=4,
							  OUTPUT_VOL=8,
							  CON_NODES=16,
							  CON_PATH=32,
  							  KEY_FRAME = 64,
							  CON_TRAJECTORY = 128};
  /**
   * @class AniEditDMRender render the data of class AniEditDM.
   * 
   */
  class AniEditDMRender: public SelfRenderEle{
	
  public:
	AniEditDMRender(pAniEditDM data_model, int render_type, const string title);
	void render();
	void draw()const;
	void updateConShpereRadius();

	void setRenderType(const int type){
	  render_type = ( type >= 0 ? type : NOTHING );
	}
	void addRenderType(const ANIEDITDM_RENDER_TYPE type){
	  render_type = render_type | type;
	}
	void removeRenderType(const ANIEDITDM_RENDER_TYPE type){
	  render_type = render_type & (~type);
	}
	void toggleRenderType(const ANIEDITDM_RENDER_TYPE type){
	  includeRenderType(type) ? removeRenderType(type) : addRenderType(type);
	}
	bool includeRenderType(const ANIEDITDM_RENDER_TYPE type)const{
	  return (render_type&type) != NOTHING;
	}
	const int getRenderType()const{
	  return render_type;
	}

  protected:
	void updateFrameNumText();
	void drawInputObj()const;
	void drawInputVol()const;
	void drawOutputObj()const;
	void drawKeyframeObj()const;
	void drawOutputVol()const;
	void drawConNodes()const;
	void drawConPath()const;
	void drawConTraj()const;

	bool hasVolMesh()const{
	  return (data_model != NULL && data_model->getVolMesh() != NULL);
	}
	pVolumetricMesh_const getVolMesh()const{
	  return data_model->getVolMesh();
	}
	
  private:
	pAniEditDM data_model;
	LSW_BASE_UI::VolNodeGroupRender node_group_render;
	int render_type;
	float input_obj_fAmbient[4];
	float input_obj_fDiffuse[4];
	string title;
  };
  
  typedef boost::shared_ptr<AniEditDMRender> pAniEditDMRender;
  
}//end of namespace

#endif /*_ANIEDITDMRENDER_H_*/
