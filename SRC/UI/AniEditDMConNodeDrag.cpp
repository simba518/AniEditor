#include <GL/glew.h>
#include <GL/gl.h>
#include "DragGroupSelDraw.h"
#include "AniEditDMConNodeDrag.h"
using namespace LSW_ANI_EDIT_UI;

AniEditDMConNodeDrag::AniEditDMConNodeDrag(pQGLViewerExt viewer, pAniEditDM dm):
  viewer(viewer),data_model(dm){
  
  drag_point.setZero();
  draw_mouse_circle = false;
  draw_drag_circle = false;
  mouse_point << -1e6,-1e6,-1e6;
  draged_point << -1e6,-1e6,-1e6;

  const double dragged_color[4] = {0.0f, 0.0f, 0.6f, 1.0f};
  con_node_render.setColor(dragged_color);
}

void AniEditDMConNodeDrag::drawWithNames ()const{

  pTetMesh_const vol_mesh = data_model->getVolMesh();
  if( vol_mesh != NULL && data_model->currentFrameNum() >= 0){
	
	const vector<set<int> > con_groups = data_model->getConNodes();
	const VectorXd &vol_u = data_model->getVolFullU();
	const double radius = data_model->getMaxRadius()/78.0f;
	assert_gt (radius,0.0f);
	if (vol_u.size() > 0){
	  LSW_SIM_UI::DragGroupSelDraw::drawAllGroupsWithShpere(vol_mesh,con_groups,&vol_u[0], radius);
	}
  }
}

void AniEditDMConNodeDrag::render(){
  
  const double radius = data_model->getMaxRadius()/78.0f;
  if (radius > 0.0f){
	  con_node_render.setSphereRadius(radius);
	}
  this->draw();
}

void AniEditDMConNodeDrag::draw()const{
  
  if (data_model != NULL && data_model->currentFrameNum() >=0 
	  && dragger.hasDragedGroup()){

	const set<int> &dragged_nodes = dragger.getSelGroupNods();
	pTetMesh_const vol_mesh = data_model->getVolMesh();
	const VectorXd &vol_u = data_model->getVolFullU();
	
	if (vol_u.size() >= 3){
	  con_node_render.draw(vol_mesh, dragged_nodes, &vol_u[0], DRAW_SHPERE);
	}
  }

  if (draw_mouse_circle){
	drawMouseCircle();
  }

  if (draw_drag_circle){
	drawDragPointCircle();
  }
}

void AniEditDMConNodeDrag::getDragedPoint(double point[3])const{
  point[0] = drag_point[0];
  point[1] = drag_point[1];
  point[2] = drag_point[2];
}

void AniEditDMConNodeDrag::selectDragEle(int sel_group_id){

  dragger.setDragGroup(sel_group_id);
  if (data_model){

	const vector<set<int> > groups = data_model->getConNodes();
	pTetMesh_const vol_mesh = data_model->getVolMesh();
	if( vol_mesh && sel_group_id < (int)groups.size() ) {

	  drag_point = barycentersOfGroup(vol_mesh->nodes(),groups[sel_group_id]);
	  const VectorXd &u = data_model->getVolFullU();
	  drag_point += getBaryCenter(groups[sel_group_id],u);
	}
  }
  updateConPos();
}

const Vector3d AniEditDMConNodeDrag::getBaryCenter
(const set<int> &one_group, const VectorXd &u)const{

  Vector3d bary_p;
  bary_p[0] = 0;
  bary_p[1] = 0;
  bary_p[2] = 0;
  
  std::set<int>::const_iterator it = one_group.begin();
  for (; it != one_group.end(); it++ ){

	Vector3d v;
	v[0] = u[(*it)*3+0];
	v[1] = u[(*it)*3+1];
	v[2] = u[(*it)*3+2];
	bary_p[0] += v[0];
	bary_p[1] += v[1];
	bary_p[2] += v[2];
  }

  if (one_group.size() > 0){

	bary_p[0] /= one_group.size();
	bary_p[1] /= one_group.size();
	bary_p[2] /= one_group.size();
  }

  return bary_p;
}

void AniEditDMConNodeDrag::drawMouseCircle()const{

  if (viewer){
	const double color[3] = {1,0,0};
	const double wx = mouse_point[0];
	const double wy = mouse_point[1];
	const double wz = mouse_point[2];
	const Vec v = viewer->getScreenCoords(wx,wy,wz);
	draw2DCircle(viewer->width(), viewer->height(), v[0], v[1], color);
  }
}

void AniEditDMConNodeDrag::drawDragPointCircle()const{

  if (viewer){

	const double color[3] = {0,0,1};
	const double wx = draged_point[0];
	const double wy = draged_point[1];
	const double wz = draged_point[2];
	const Vec v = viewer->getScreenCoords(wx,wy,wz);
	draw2DCircle(viewer->width(), viewer->height(), v[0], v[1], color);
  }
}

void AniEditDMConNodeDrag::draw2DCircle(const double screen_width,
										const double screen_height,
										const double screen_x, 
										const double screen_y, 
										const double color[3])const{
  
  assert_gt (screen_width, 0);
  assert_gt (screen_height, 0);

  glDisable(GL_LIGHTING);

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0.0f,screen_width,0.0f,screen_height);

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  
  glColor3f(color[0],color[1],color[2]);

  glTranslated(screen_x, screen_height-screen_y,0.0f);
  gluDisk(gluNewQuadric(), 18,20, 100,100);

  glPopMatrix();
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);

  glEnable(GL_LIGHTING);
}

void AniEditDMConNodeDrag::updateConPos(){

  if (data_model != NULL){
	data_model->updateConPos( dragger.getUcVec() );
  }

  if (viewer && data_model->getVolMesh() && dragger.hasDragedGroup()){

	const set<int> &dragged_nodes = dragger.getSelGroupNods();
	pTetMesh_const vol_mesh = data_model->getVolMesh();
	const Vector3d bary_p = barycentersOfGroup(vol_mesh->nodes(),dragged_nodes);

	// initialize the center of the circle of the draged point.
	const VectorXd &vol_u = data_model->getVolFullU();
	const Vector3d cu = getBaryCenter(dragged_nodes, vol_u);
	const double w_x = bary_p[0]+cu[0];
	const double w_y = bary_p[1]+cu[1];
	const double w_z = bary_p[2]+cu[2];
	draged_point << w_x, w_y, w_z;

	// initialize the center of the circle of the mouse.
	double p[3];
	getDragedPoint(p);
	const Vec screen_pos = viewer->getScreenCoords(p[0],p[1],p[2]);
	const double z_deepth = screen_pos[2];
	const Vec world_pos = viewer->getWorldCoords(viewer->getMousePos(),z_deepth);
	mouse_point << world_pos[0],world_pos[1],world_pos[2];
  }

  if (viewer != NULL){
	viewer->update();
  }
}

Vector3d AniEditDMConNodeDrag::barycentersOfGroup(const VVec3d &nodes, const set<int> &g)const{

  Vector3d bary_p;
  bary_p.setZero();
  set<int>::const_iterator it = g.begin();
  for (; it != g.end(); ++it){
	assert_in(*it,0,nodes.size()-1);
	bary_p += nodes[*it];
  }
  if (g.size() > 0)
	bary_p /= g.size();
  return bary_p;
}
