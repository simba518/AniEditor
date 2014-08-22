#ifndef _MANIPULATEOP_H_
#define _MANIPULATEOP_H_

#include <QGLViewerExt.h>
#include <SelectCtrl.h>
#include <Manipulatoion.h>
#include <Selectable.h>
#include <DragGroupSelDraw.h>
#include "AniEditDM.h"
using namespace Eigen;
using namespace UTILITY;
using namespace QGLVEXT;

namespace ANI_EDIT_UI{
  
  /**
   * @class ManipulateOP
   * 
   */
  class ManipulateOP:public LocalframeManipulatoionExt{

	Q_OBJECT
	
  public:
	ManipulateOP(pQGLViewerExt v,pAniEditDM dm):
	  LocalframeManipulatoionExt(v),data_model(dm){
	  dragged_group_id = -1;
	}

	// select
	int totalEleNum()const{
	  return data_model != NULL ? data_model->getConNodes().size() : 0;
	}
	void drawWithNames()const{
	  pTetMesh_const vol_mesh = data_model->getVolMesh();
	  if( vol_mesh != NULL ){
		const vector<set<int> > con_groups = data_model->getConNodes();
		const VectorXd &vol_u = data_model->getVolFullU();
		double x=0,y=0,z=0;
		if (data_model){
		  data_model->getSceneTranslate(x,y,z);
		  glTranslated(x,y,z);
		}
		DragGroupSelDraw::drawAllGroupsWithPoints(vol_mesh,con_groups,&vol_u[0]);
		glTranslated(-x,-y,-z);

	  }
	}
	void select(const vector<int> &sel_ids){

	  if(sel_ids.size()>0 && data_model && data_model->getVolMesh()){
		assert_in(sel_ids[0],0,totalEleNum());
		dragged_group_id = sel_ids[0];
		manipulated_pos = data_model->getXc(sel_ids[0]);
		assert_gt(manipulated_pos.size(),0);
	  }
	}

	Matrix<double,3,-1> &getCurrentPositions(){
	  return manipulated_pos;
	}
	const Matrix<double,3,-1> &getCurrentPositions()const{
	  return manipulated_pos;
	}
	void applyTransform(){
	  LocalframeManipulatoionExt::applyTransform();
	  if(data_model && dragged_group_id >= 0){
		data_model->updateXc(manipulated_pos,dragged_group_id);
	  }
	}

  public slots:
	void updateCameraTrans(){
	  double x=0,y=0,z=0;
	  if (data_model)
		data_model->getSceneTranslate(x,y,z);
	  cout << "updateCameraTrans: "<< x << ", "<< y << ", " << z << endl;
	  this->translateCamera(x,y,z);
	}

  private:
	pAniEditDM data_model;
	Matrix<double,3,-1> manipulated_pos;
	int dragged_group_id;
  };
  typedef boost::shared_ptr<ManipulateOP> pManipulateOP;
  
}//end of namespace

#endif /*_MANIPULATEOP_H_*/
