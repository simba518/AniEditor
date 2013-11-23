#ifndef _ANIEDITDMRENDER_H_
#define _ANIEDITDMRENDER_H_

#include <boost/shared_ptr.hpp>
#include <SelfRenderEle.h>
#include <TextForRender.h>
#include <VolNodeGroupRender.h>
#include <QObject>
#include <QGLViewerExt.h>
#include "AniEditDM.h"
#include "AniEditDMRender.h"
using namespace QGLVEXT;

namespace ANI_EDIT_UI{

  enum ANIEDITDM_RENDER_TYPE {NOTHING=0,
							  INPUT_OBJ=1,
							  INPUT_VOL=2,
							  OUTPUT_OBJ=4,
							  OUTPUT_VOL=8,
							  CON_NODES=16,
							  CON_PATH=32,
							  CON_TRAJECTORY = 128};


  /// @class AniEditDMRender render the data of class AniEditDM.
  class AniEditDMRender: public SelfRenderEle, public TextForRender{
	
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
	void drawOutputVol()const;
	void drawConNodes()const;

	bool hasVolMesh()const{
	  return (data_model != NULL && data_model->getVolMesh() != NULL);
	}
	pTetMesh_const getVolMesh()const{
	  return data_model->getVolMesh();
	}
	
  private:
	pAniEditDM data_model;
	UTILITY::VolNodeGroupRender node_group_render;
	int render_type;
	ObjMtl input_obj_mtl;
	string title;
  };
  typedef boost::shared_ptr<AniEditDMRender> pAniEditDMRender;


  /// @class AniEditDMRenderCtrl control the rendering of AniEditDM class.
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

#endif /*_ANIEDITDMRENDER_H_*/
