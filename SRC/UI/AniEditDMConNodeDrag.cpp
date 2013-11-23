#include <GL/glew.h>
#include <GL/gl.h>
#include <DragGroupSelDraw.h>
#include "AniEditDMConNodeDrag.h"
using namespace ANI_EDIT_UI;
using namespace UTILITY;

AniEditDMConNodeDrag::AniEditDMConNodeDrag(pQGLViewerExt viewer, pAniEditDM dm):
  viewer(viewer),data_model(dm){
  draggedGroupId = -1;
  dragged_point_start.setZero();
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

  point[0] = dragged_point_start[0];
  point[1] = dragged_point_start[1];
  point[2] = dragged_point_start[2];
}

void AniEditDMConNodeDrag::selectDragEle(int sel_group_id){

  if(data_model){
	assert_in(sel_group_id,0,data_model->getConNodes().size()-1);
	draggedGroupId = sel_group_id;
	const set<int> groups = data_model->getConNodes()[sel_group_id];
	pTetMesh_const restVolMesh = data_model->getVolMesh();
	assert_gt(groups.size(),0);
	const int nodeI = *groups.begin();
	dragged_point_start = restVolMesh->nodes()[nodeI]+data_model->getUc(sel_group_id).col(0);
  }
}

void AniEditDMConNodeDrag::dragTo(double x,double y,double z){
  if(data_model){
	Vector3d u;
	u << x,y,z;
	u -= dragged_point_start;
	data_model->updateConPos(u,draggedGroupId);
  }
}
