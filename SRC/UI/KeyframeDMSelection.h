#ifndef _KEYFRAMEDMSELECTION_H_
#define _KEYFRAMEDMSELECTION_H_

#include <QObject>
#include <SelectCtrl.h>
#include <Selectable.h>
#include <QGLViewerExt.h>
#include "KeyframeDM.h"
using namespace QGLVEXT;
using namespace UTILITY;

namespace ANI_EDIT_UI{
  
  /**
   * @class KeyframeDMSelection selection constraint nodes for keyframe data model.
   * 
   */
  class KeyframeDMSelection:public Selectable, public SelectObserver{
	
  public:
	KeyframeDMSelection(pKeyframeDM dm):data_model(dm){}
	int totalEleNum ()const{
	  return hasVolMesh() ? data_model->getVolMesh()->nodes().size():0;
	}
	void drawWithNames ()const{drawVertice();}
	void addSelection(const vector<int> &sel_ids){
	  assert(data_model->isKeyframe()>=0);
	  if (data_model != NULL) data_model->addConNodes(sel_ids);
	}
	void removeSelection(const vector<int> &sel_ids){
	  assert(data_model->isKeyframe()>=0);
	  if (data_model != NULL) data_model->rmConNodes(sel_ids);
	}
	
  protected:
	bool hasVolMesh()const{
	  return (data_model != NULL) && (data_model->getVolMesh() != NULL&&data_model->isKeyframe()>=0);
	}
	void drawVertice()const{

	  if( this->hasVolMesh() && data_model->currentFrameNum()>=0){
		pTetMesh_const vol_mesh = data_model->getVolMesh();
		const VectorXd &vol_u = data_model->getKeyU();
		glFlush();
		for (int i=0; i<vol_mesh->nodes().size(); i++){
		  const Vector3d &v = vol_mesh->nodes()[i];
		  glPushName(i);
		  glBegin(GL_POINTS);
		  glVertex3d(v[0]+vol_u[i*3+0], v[1]+vol_u[i*3+1], v[2]+vol_u[i*3+2]);
		  glEnd();
		  glPopName();
		}
	  }
	}
	
  private:
	pKeyframeDM data_model;
  };
  typedef boost::shared_ptr<KeyframeDMSelection> pKeyframeDMSelection;

  class KeyframeDMSelectionCtrl: public QObject{

	Q_OBJECT
	
  public:
	KeyframeDMSelectionCtrl(pQGLViewerExt viewer, pKeyframeDM data_model){
	  vol_vert_sel=pKeyframeDMSelection(new KeyframeDMSelection(data_model));
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
	pKeyframeDMSelection vol_vert_sel;
  };
  typedef boost::shared_ptr<KeyframeDMSelectionCtrl> pKeyframeDMSelectionCtrl;
  
}//end of namespace

#endif /*_KEYFRAMEDMSELECTION_H_*/
