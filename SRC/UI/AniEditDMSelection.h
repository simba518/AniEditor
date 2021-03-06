#ifndef _ANIEDITDMSELECTION_H_
#define _ANIEDITDMSELECTION_H_

#include <boost/shared_ptr.hpp>
#include <QObject>
#include <SelectCtrl.h>
#include <Selectable.h>
#include "AniEditDM.h"
using namespace QGLVEXT;

namespace ANI_EDIT_UI{

  /**
   * @class AniEditDMVolVertexSel implement the drawWithNames() function to
   * select the vertex of the volume mesh of current frame of AniEditDM data
   * model.
   */
  class AniEditDMVolVertexSel:public Selectable, public SelectObserver{
	
  public: 
	AniEditDMVolVertexSel(pAniEditDM dm):data_model(dm){}
	int totalEleNum ()const{
	  return hasVolMesh() ? data_model->getVolMesh()->nodes().size():0;
	}
	void drawWithNames ()const{drawVertice();}
	void addSelection(const vector<int> &sel_ids){
	  if (sel_ids.size() > 0){
		// vector<int> first_node;
		// first_node.push_back(sel_ids[0]);
		// if (data_model != NULL) data_model->addConNodes(first_node);
		if (data_model != NULL) data_model->addConNodes(sel_ids);
	  }
	}
	void removeSelection(const vector<int> &sel_ids){
	  if (data_model != NULL) data_model->rmConNodes(sel_ids);
	}

  protected:
	bool hasVolMesh()const{
	  return (data_model != NULL) && (data_model->getVolMesh() != NULL);
	}
	void drawVertice()const{
	  
	  if( this->hasVolMesh() && data_model->currentFrameNum()>=0){

		pTetMesh_const vol_mesh = data_model->getVolMesh();
		const VectorXd &vol_u = data_model->getVolFullU();

		double x=0,y=0,z=0;
		if (data_model){
		  data_model->getSceneTranslate(x,y,z);
		  glTranslated(x,y,z);
		}
		glFlush();
		for (int i=0; i<vol_mesh->nodes().size(); i++){
		  const Vector3d &v = vol_mesh->nodes()[i];
		  glPushName(i);
		  glBegin(GL_POINTS);
		  glVertex3d(v[0]+vol_u[i*3+0], v[1]+vol_u[i*3+1], v[2]+vol_u[i*3+2]);
		  glEnd();
		  glPopName();
		}
		glTranslated(-x,-y,-z);

	  }
	}
	
  private:
	pAniEditDM data_model;
  };
  typedef boost::shared_ptr<AniEditDMVolVertexSel> pAniEditDMVolVertexSel;

  
  /// @class AniEditDMSelCtrl control the select operation of AniEditDM object.
  class AniEditDMSelCtrl: public QObject{

	Q_OBJECT
	
  public:
	AniEditDMSelCtrl(pQGLViewerExt viewer, pAniEditDM data_model){
	  vol_vert_sel=pAniEditDMVolVertexSel(new AniEditDMVolVertexSel(data_model));
	  sel_ctrl = pSelectCtrl(new SelectCtrl(viewer,vol_vert_sel,10));
	  sel_ctrl->setObserver(vol_vert_sel);
	}
	
  public slots:
	void togglePrintSelEle(){ sel_ctrl->togglePrintSelEle(); }
	void enableSelOp(){	sel_ctrl->enableSelOp(); }
	void disableSelOp(){ sel_ctrl->disableSelOp(); }
	bool toggleSelOp(){	return sel_ctrl->toggleSelOp();	}

  private:
	pSelectCtrl sel_ctrl;
	pAniEditDMVolVertexSel vol_vert_sel;
  };
  typedef boost::shared_ptr<AniEditDMSelCtrl> pAniEditDMSelCtrl;
  
}//end of namespace

#endif /* _ANIEDITDMSELECTION_H_ */
