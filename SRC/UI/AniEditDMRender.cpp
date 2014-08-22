#include <WindowsHeaders.h>
#include <GL/gl.h>
#include <boost/lexical_cast.hpp>
#include "AniEditDMRender.h"
#include <MeshRender.h>
using namespace ANI_EDIT_UI;
using namespace UTILITY;

AniEditDMRender::AniEditDMRender(pAniEditDM dm, pQGLViewerExt viewer,int type, const string title):
  data_model(dm),_viewer(viewer),render_type(type),title(title){

  input_obj_mtl.setDefault();
  add_obj_mtl.setDefault();
  add_obj_mtl.ambient[0] = 0.6f;
  add_obj_mtl.diffuse[0] = 0.6f;
  const double con_group_color[4] = {0.6f, 0.0f, 0.0f, 1.0f};
  node_group_render.setColor(con_group_color);

  tran.setZero();
  rot_axi << 1.0f,0.0f,0.0f;
  angle = 0.0f;
}

void AniEditDMRender::render(){

  updateFrameNumText();
  updateConShpereRadius();
  this->draw();
}

void AniEditDMRender::draw()const{

  double x=0,y=0,z=0;
  if (data_model){
	data_model->getSceneTranslate(x,y,z);
	glTranslated(x,y,z);
  }

  if (render_type & OUTPUT_OBJ)
  	drawOutputObj();
  if (render_type & OUTPUT_VOL)
  	drawOutputVol();
  if (render_type & INPUT_OBJ)
  	drawInputObj();
  if (render_type & INPUT_VOL)
  	drawInputVol();
  if (render_type & CON_NODES)
  	drawConNodes();
  if (render_type & CON_PATH)
  	drawConPath();
  if (render_type & ADDITIONAL_ANI)
	drawAdditionalAniObj();
  if (render_type & SCENE)
	drawScene();

  glTranslated(-x,-y,-z);

  if (render_type & GROUND)
	drawGround();

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

  	// if(data_model->getInterpolator()){
  	//   const string interp_method = data_model->getInterpolator()->getName();
  	//   this->update(interp_method,x,y+90);
  	//   const string warp = data_model->getInterpolator()->isUseWarp()? string("warp"):string("unwarp");
  	//   this->update(warp,x,y+90+45);
  	// }
  }
}

void AniEditDMRender::drawInputObj()const{
  if (data_model != NULL){
	// const VectorXd &vol_u = data_model->getInputU();
	// cout<< "input: " << (vol_u.segment(1555*3,3)+data_model->getVolMesh()->nodes()[1555]).transpose() << endl;
	UTILITY::draw(data_model->getInputObjMesh(),input_obj_mtl);
  }
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

void AniEditDMRender::drawAdditionalAniObj()const{
  if (data_model != NULL && data_model->hasAdditionalAni())
	UTILITY::draw(data_model->getAdditionalAniObjMesh(),add_obj_mtl);
}

void AniEditDMRender::drawConPath()const{

  if ( hasVolMesh() ){

	const VVec3d& rest_u = getVolMesh()->nodes();
	const set<pPartialConstraints>& parcon = data_model->getAllConNodes().getPartialConSet();

	glBegin(GL_POINTS);
	BOOST_FOREACH(pPartialConstraints_const con, parcon){

	  if (con){
		const vector<set<int> > &group = con->getConNodesSet();
		const Matrix<double,3,-1> &con_u = con->getPc();
		int p = 0;
		for (int i = 0; i < group.size(); ++i){

		  BOOST_FOREACH(const int &node_i, group[i]){
			assert_in(p,0,con_u.cols()-1);
			assert_in(node_i,0,rest_u.size()-1);
			const Vector3d v = con_u.col(p)+rest_u[node_i];
			glVertex3d(v[0],v[1],v[2]);
			p ++;
		  }
		}
	  }
	}
	glEnd();  
  }

  {// @todo to remove, save constraint nodes.
	// const string d = "/home/simba/Workspace/AnimationEditor/Data/flower_box_casa/";
	// const string ball_f = d+"model/ball.obj";
	// Objmesh ball;
	// bool succ = ball.load(ball_f);
	// ERROR_LOG_COND("failed to load " << ball_f,succ);
	
	// const VVec3d& rest_x = getVolMesh()->nodes();
	// const set<pPartialConstraints>& parcon = data_model->getAllConNodes().getPartialConSet();
	// const Eigen::VectorXd vertices = ball.getVerts();
	// const Eigen::Vector3d cen = ball.getCenter();
	// BOOST_FOREACH(pPartialConstraints_const con, parcon){
	//   if (con){
	// 	const vector<set<int> > &group = con->getConNodesSet();
	// 	const Matrix<double,3,-1> &con_u = con->getPc();
	// 	int p = 0;
	// 	const int f = con->getFrameId();
	// 	for (int i = 0; i < group.size(); ++i){
	// 	  BOOST_FOREACH(const int &node_i, group[i]){
	// 		// const Vector3d pp = rest_x[node_i]+data_model->getVolFullU(f).segment<3>(node_i*3)-cen;
	// 		const Vector3d pp = rest_x[node_i]+con_u.col(p)-cen;
	// 		p ++;
	// 		VectorXd pos = vertices;
	// 		for (int i = 0; i < pos.size()/3; i++)
	// 		  pos.segment<3>(i*3) += pp;
	// 		ball.setVerts(pos);
	// 		const string save_to = "./tempt/con_node_"+TOSTR(node_i)+"_frame_"+TOSTR(f)+"_ball.obj";
	// 		cout << "save to : " << save_to << endl;
	// 		succ = ball.write(save_to);
	// 		ERROR_LOG_COND("failed to write: " << save_to, succ);
	// 	  }
	// 	}
	//   }
	// }
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

void AniEditDMRender::updateConShpereRadius(){
  
  if (data_model){
	const double radius = data_model->getMaxRadius()/50.0f;
	if (radius > 0.0f){
	  node_group_render.setSphereRadius(radius);
	}
  }
}

void AniEditDMRender::drawGround()const{

  if (_viewer){
	glColor3f(0.1,0.1,0.1);
	glLineWidth(1.0f);
	glPushMatrix();
	glTranslatef(tran[0], tran[1], tran[2]);
	glRotatef(angle, rot_axi[0], rot_axi[1], rot_axi[2]);
	QGLViewer::drawGrid(10,300);
	glPopMatrix();
  }
}

void AniEditDMRender::drawScene()const{
  
  if (data_model != NULL && data_model->getScene())
	UTILITY::draw(data_model->getScene());
}
