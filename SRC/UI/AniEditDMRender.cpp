#include <GL/gl.h>
#include <boost/lexical_cast.hpp>
#include <ConNodeSetDrawer.h>
#include "AniEditDMRender.h"
#include <MeshRender.h>
using namespace ANI_EDIT_UI;
using namespace UTILITY;

AniEditDMRender::AniEditDMRender(pAniEditDM dm,int type, const string title)
  :data_model(dm),render_type(type),title(title){

  input_obj_mtl.setDefault();
  const double con_group_color[4] = {0.6f, 0.0f, 0.0f, 1.0f};
  node_group_render.setColor(con_group_color);
}

void AniEditDMRender::render(){

  updateFrameNumText();
  updateConShpereRadius();
  this->draw();
}

void AniEditDMRender::draw()const{

  if (render_type & OUTPUT_OBJ){
  	drawOutputObj();
  }
  if (render_type & OUTPUT_VOL){
  	drawOutputVol();
  } 
  if (render_type & INPUT_OBJ){
  	drawInputObj();
  }
  if (render_type & INPUT_VOL){
  	drawInputVol();
  }
  if (render_type & CON_NODES){
  	drawConNodes();
  }
  if (render_type & CON_PATH){
  	drawConPath();
  }
  if (render_type & CON_TRAJECTORY){
  	drawConTraj();
  }
}

void AniEditDMRender::updateFrameNumText(){

  if (data_model != NULL){
	
  	const int x = 60;
  	const int y = 70;
  	this->update(this->title,x,y);
	
  	const int f = data_model->currentFrameNum();
  	const string frame_num = boost::lexical_cast<string>(f);
  	const string T = boost::lexical_cast<string>( data_model->totalFrameNum()-1);
  	const string text = string("Timestep  ") + frame_num + string(" / ")+T;
  	this->update(text,x,y+45);

  	if(data_model->getInterpolator()){
  	  const string interp_method = data_model->getInterpolator()->getName();
  	  this->update(interp_method,x,y+90);
  	  const string warp = data_model->getInterpolator()->isUseWarp()? string("warp"):string("unwarp");
  	  this->update(warp,x,y+90+45);
  	}
  }
}

void AniEditDMRender::drawInputObj()const{
  if (data_model != NULL)
	UTILITY::draw(data_model->getInputObjMesh(),input_obj_mtl);
}

void AniEditDMRender::drawInputVol()const{

  if (data_model != NULL && data_model->getVolMesh() != NULL){
	pTetMesh_const vol_mesh = data_model->getVolMesh();
	const double *u = NULL;
	if (data_model->currentFrameNum() >=0){
	  const VectorXd &vol_u = data_model->getInputU();
	  if(vol_u.size() > 0) u = &(vol_u[0]);
	}
	UTILITY::draw(vol_mesh,u);
  }
}

void AniEditDMRender::drawOutputObj()const{
  if (data_model != NULL)
	UTILITY::draw(data_model->getOutputObjMesh());
}

void AniEditDMRender::drawOutputVol()const{
  
  if (data_model != NULL && data_model->getVolMesh() != NULL){
	pTetMesh_const vol_mesh = data_model->getVolMesh();
	if (data_model->currentFrameNum() >=0){
	  const VectorXd &vol_u = data_model->getVolFullU();
	  if(vol_u.size() > 0){
		assert_eq(vol_u.size(),vol_mesh->nodes().size()*3);
		UTILITY::draw(vol_mesh,&(vol_u[0]));
	  }
	}else{
	  UTILITY::draw(vol_mesh);
	}
  }
}

void AniEditDMRender::drawConNodes()const{

  if ( hasVolMesh() && data_model->currentFrameNum() >= 0){

	pTetMesh_const vol_mesh = getVolMesh();
	const vector<set<int> > groups = data_model->getConNodes();
	const VectorXd &vol_u = data_model->getVolFullU();
	if (vol_u.size() >= 3){
	  node_group_render.draw(vol_mesh, groups,&vol_u[0],UTILITY::DRAW_POINT);
	}
  }
}

void AniEditDMRender::drawConPath()const{

  if (data_model != NULL){
	const ConNodesOfFrameSet &c_nodes=data_model->getAllConNodes();
	const double radius = data_model->getMaxRadius()/50.0f;
	ConNodeSetDrawer::draw(c_nodes.getConNodeGroups(),radius);
  }
}

void AniEditDMRender::updateConShpereRadius(){
  
  if (data_model){
	const double radius = data_model->getMaxRadius()/50.0f;
	if (radius > 0.0f){
	  node_group_render.setSphereRadius(radius);
	}
  }
}

void AniEditDMRender::drawConTraj()const{
  
  if(data_model){
	const vector<VectorXd> &con_trajs = data_model->getConNodeTraj();
	
	for (size_t i = 0; i < con_trajs.size(); ++i){
	  const VectorXd &p = con_trajs[i];
	  glDisable(GL_LIGHTING);
	  glColor3d(1,0,0);
	  glLineWidth(1.0f);
	  glBegin(GL_LINE_STRIP);
	  for (int i = 0; i < p.size()/3; ++i){
		const double x = p[i*3+0];
		const double y = p[i*3+1];
		const double z = p[i*3+2];
		glVertex3f(x,y,z);
	  }
	  glEnd();
	  glEnable(GL_LIGHTING);
	}
  }
}
