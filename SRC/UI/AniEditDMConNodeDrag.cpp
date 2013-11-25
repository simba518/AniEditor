#include <WindowsHeaders.h>
#include <GL/glew.h>
#include <GL/gl.h>
#include <DragGroupSelDraw.h>
#include "AniEditDMConNodeDrag.h"
using namespace ANI_EDIT_UI;
using namespace UTILITY;

AniEditDMConNodeDrag::AniEditDMConNodeDrag(pQGLViewerExt viewer, pAniEditDM dm):
  viewer(viewer),data_model(dm){
  draggedGroupId = -1;
  initial_dragged_point.setZero();
  initial_displacement.setZero();
  const double dragged_color[4] = {0.0f, 0.0f, 0.6f, 1.0f};
  con_node_render.setColor(dragged_color);
}

void AniEditDMConNodeDrag::drawWithNames ()const{

  pTetMesh_const vol_mesh = data_model->getVolMesh();
  if( vol_mesh != NULL && data_model->currentFrameNum() >= 0){
	
  	const vector<set<int> > con_groups = data_model->getConNodes();
  	const VectorXd &vol_u = data_model->getVolFullU();
  	const double radius = data_model->getMaxRadius()/78.0f;
  	if (vol_u.size() > 0 && radius > 0.0f){
  	  DragGroupSelDraw::drawAllGroupsWithShpere(vol_mesh,con_groups,&vol_u[0], radius);
  	}
  }
}

void AniEditDMConNodeDrag::render(){
  
  const double radius = data_model->getMaxRadius()/78.0f;
  if (radius > 0.0f)
	con_node_render.setSphereRadius(radius);
  this->draw();
}

void AniEditDMConNodeDrag::draw()const{
  
  if (data_model!=NULL&&data_model->currentFrameNum()>=0&&hasDraggedGroup()){
  	pTetMesh_const vol_mesh = data_model->getVolMesh();
	const set<int> dragged_nodes = data_model->getConNodes()[draggedGroupId];
  	const VectorXd &vol_u = data_model->getVolFullU();
  	if (vol_u.size() >= 3)
  	  con_node_render.draw(vol_mesh, dragged_nodes, &vol_u[0], DRAW_SHPERE);
  }
}

void AniEditDMConNodeDrag::getDragedPoint(double point[3])const{

  point[0] = initial_dragged_point.col(0)[0];
  point[1] = initial_dragged_point.col(0)[1];
  point[2] = initial_dragged_point.col(0)[2];
}

void AniEditDMConNodeDrag::selectDragEle(int sel_group_id){

  if(data_model){
	assert_in(sel_group_id,0,data_model->getConNodes().size()-1);
	draggedGroupId = sel_group_id;
	const set<int> groups = data_model->getConNodes()[sel_group_id];
	pTetMesh_const restVolMesh = data_model->getVolMesh();
	assert_gt(groups.size(),0);

	initial_displacement = data_model->getUc(sel_group_id);
	assert_eq(initial_displacement.cols(),groups.size());
	initial_dragged_point = initial_displacement.col(0);
	initial_dragged_point += restVolMesh->nodes()[*groups.begin()];
  }
}

void AniEditDMConNodeDrag::dragTo(double x,double y,double z){

  if(data_model){
	double disp[3];
	disp[0] = x-initial_dragged_point[0];
	disp[1] = y-initial_dragged_point[1];
	disp[2] = z-initial_dragged_point[2];

	Matrix<double,3,-1> u = initial_displacement;
	for (int i = 0; i < u.cols(); ++i){
	  u.col(i)[0] += disp[0];
	  u.col(i)[1] += disp[1];
	  u.col(i)[2] += disp[2];
	}
	data_model->updateConPos(u,draggedGroupId);
  }
  if(viewer)
	viewer->update();
}
