#ifndef _ANIEDITDMRENDER_H_
#define _ANIEDITDMRENDER_H_

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
							  CON_TRAJECTORY = 128,
							  ADDITIONAL_ANI = 256,
							  GROUND = 512,
							  SCENE = 1024};

  /// @class AniEditDMRender render the data of class AniEditDM.
  class AniEditDMRender: public SelfRenderEle, public TextForRender{
	
  public:
	AniEditDMRender(pAniEditDM data_model, pQGLViewerExt viewer, int render_type, const string title);
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
	void setGround(const Vector3f &tran, const Vector3f &rot_axi, const double &angle){
	  this->tran = tran;
	  this->rot_axi = rot_axi;
	  this->angle = angle;
	}

  protected:
	void updateFrameNumText();
	void drawInputObj()const;
	void drawInputVol()const;
	void drawOutputObj()const;
	void drawOutputVol()const;
	void drawAdditionalAniObj()const;
	void drawConNodes()const;
	void drawConPath()const;
	void drawGround()const;
	void drawScene()const;

	bool hasVolMesh()const{
	  return (data_model != NULL && data_model->getVolMesh() != NULL);
	}
	pTetMesh_const getVolMesh()const{
	  return data_model->getVolMesh();
	}
	
  private:
	pAniEditDM data_model;
	pQGLViewerExt _viewer;
	UTILITY::VolNodeGroupRender node_group_render;
	int render_type;
	ObjMtl input_obj_mtl;
	ObjMtl add_obj_mtl;
	string title;
	// ground
	Vector3f tran, rot_axi;
	double angle;
  };
  typedef boost::shared_ptr<AniEditDMRender> pAniEditDMRender;

  /// @class AniEditDMRenderCtrl control the rendering of AniEditDM class.
  class AniEditDMRenderCtrl: public QObject{
	
	Q_OBJECT
	
	public:
	AniEditDMRenderCtrl(pQGLViewerExt viewer,pAniEditDM dm, const string title = "Editing"):viewer(viewer){
  
	  model_render=pAniEditDMRender(new AniEditDMRender(dm,viewer,OUTPUT_OBJ|CON_NODES|CON_TRAJECTORY|ADDITIONAL_ANI|SCENE,title));
	  if(viewer != NULL){
		viewer->addSelfRenderEle(model_render);
		viewer->addTextForRender(model_render);
	  }
	}
	void setGround(const Vector3f &tran, const Vector3f &rot_axi, const double &angle){
	  model_render->setGround(tran, rot_axi, angle);
	}

  public slots:
	void toggleShowInputObj(){ toggleShowType(INPUT_OBJ); }
	void toggleShowInputVol(){ toggleShowType(INPUT_VOL); }
	void toggleShowOutputObj(){ toggleShowType(OUTPUT_OBJ); }
	void toggleShowAdditionalAni(){ toggleShowType(ADDITIONAL_ANI); }
	void toggleShowOutputVol(){ toggleShowType(OUTPUT_VOL); }
	void toggleShowConNodes(){ toggleShowType(CON_NODES); }
	void toggleShowConPath(){ toggleShowType(CON_PATH); }
	void toggleShowGround(){ toggleShowType(GROUND); }

	void showInputSeq(bool show){ showType(show, INPUT_OBJ); }
	void showConNodes(bool show){ showType(show, CON_NODES); }
	void showConPath(bool show){ showType(show, CON_PATH); }
	void showGround(const bool show){showType(show, GROUND);}

  protected:
	void showType(bool show, ANIEDITDM_RENDER_TYPE type){
	  if(show && model_render != NULL){
		model_render->addRenderType(type);
	  }else if (model_render != NULL){
		model_render->removeRenderType(type);
	  }
	  if (viewer != NULL) viewer->update();
	}
	void toggleShowType(ANIEDITDM_RENDER_TYPE type){
	  if (model_render != NULL)
		model_render->toggleRenderType(type);
	  if (viewer != NULL) viewer->update();
	}
	
  private:
	pAniEditDMRender model_render;
	pQGLViewerExt viewer;
	
  };
  typedef boost::shared_ptr<AniEditDMRenderCtrl> pAniEditDMRenderCtrl;

}//end of namespace

#endif /*_ANIEDITDMRENDER_H_*/
