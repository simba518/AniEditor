#include "AniEditDMVolVertexSel.h"
using namespace LSW_ANI_EDIT_UI;

AniEditDMVolVertexSel::AniEditDMVolVertexSel( pAniEditDM dm):
  data_model(dm){
  
}

int AniEditDMVolVertexSel::totalEleNum ()const{

  int total_ele_num = 0;
  if ( hasVolMesh() ){
	total_ele_num = getVolMesh()->nodes().size();
  }
  return total_ele_num;
}

void AniEditDMVolVertexSel::drawWithNames ()const{
  drawVertice();
}

void AniEditDMVolVertexSel::drawVertice()const{
  
  if( this->hasVolMesh() && data_model->currentFrameNum()>=0){

	pTetMesh_const vol_mesh = getVolMesh();
	const VectorXd &vol_u = data_model->getVolFullU();

	glFlush();
	for (int i=0; i<vol_mesh->nodes().size(); i++){

	  const Vector3d &v = vol_mesh->nodes()[i];
	  const double x = v[0] + vol_u[i*3+0];
	  const double y = v[1] + vol_u[i*3+1];
	  const double z = v[2] + vol_u[i*3+2];		

	  glPushName(i);
	  glBegin(GL_POINTS);
	  glVertex3d(x,y,z);
	  glEnd();
	  glPopName();
	}
  }
}
