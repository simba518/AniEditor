#include <GL/gl.h>
#include <Log.h>
#include "VolMeshDrawer.h"
using namespace LSW_BASE_SIM;

void VolMeshDrawer::render(const VolumetricMesh* volmesh,VOLMESHRENDERTYPE t,const double *u){

  if(volmesh == NULL){
	return;
  }

  switch(t){

  case POINT: drawPoint(volmesh,u);
	break;
  case EDGE: drawEdge(volmesh,u);
	break;
  case ELEMENT: drawElement(volmesh,u);
	break;
  }
}

void VolMeshDrawer::drawPoint(const VolumetricMesh* volmesh){

  if(volmesh == NULL){
	WARN_LOG("the VolumetricMesh is null");
	return;
  }

  glDisable(GL_LIGHTING);
  glPointSize(4);
  glColor3d(1, 1, 0);
  glBegin(GL_POINTS);
  for (int i=0; i<volmesh->numVertices(); i++){
	const Vec3d &v = *(volmesh->vertex(i));
	glVertex3d(v[0],v[1],v[2]);
  }
  glEnd();
  glEnable(GL_LIGHTING);
}

void VolMeshDrawer::drawPoint(const VolumetricMesh* volmesh,const double *u){
  
  if(volmesh == NULL){
	WARN_LOG("the VolumetricMesh is null");
	return;
  }

  if (u == NULL) {
	drawPoint(volmesh);
  }else{
	glDisable(GL_LIGHTING);
	glPointSize(4);
	glColor3d(1, 1, 0);
	glBegin(GL_POINTS);
	for (int i=0; i<volmesh->numVertices(); i++){
	  const Vec3d &v = *(volmesh->vertex(i));
	  double x = v[0]+u[i*3+0];
	  double y = v[1]+u[i*3+1];
	  double z = v[2]+u[i*3+2];
	  glVertex3d(x,y,z);
	}
	glEnd();
	glEnable(GL_LIGHTING);	
  }
}

void VolMeshDrawer::drawEdge(const VolumetricMesh* volmesh){
  
  ///@todo draw edges of volumetric mesh.
  
}

void VolMeshDrawer::drawEdge(const VolumetricMesh* volmesh,const double *u){
  
  if (u == NULL) {
	drawEdge(volmesh);
  }else{
	///@todo draw edges of volumetric mesh.
  }
}

void VolMeshDrawer::drawElement(const VolumetricMesh* volmesh,const double *u){

  if(volmesh == NULL){
	WARN_LOG("the VolumetricMesh is null");
	return;
  }
  
  if (u == NULL) {
	drawElement(volmesh);
  }else{

	
	glDisable(GL_LIGHTING);
	glColor3d(0.8, 0.8, 0.8);
	glBegin(GL_LINES);
	const int num_ele = volmesh->numElements();
	const int num_ele_edges = volmesh->numElementEdges();
	int edge_index[num_ele_edges*2];
	for (int el = 0; el < num_ele; ++el){

	  volmesh->getElementEdges(el,edge_index);
	  for (int i = 0; i < num_ele_edges; ++i){

		const int vi0 = edge_index[i*2];
		const int vi1 = edge_index[i*2+1];
		const Vec3d &v0 = *(volmesh->vertex(vi0));
		const Vec3d &v1 = *(volmesh->vertex(vi1));
		glVertex3d(v0[0]+u[vi0*3],v0[1]+u[vi0*3+1],v0[2]+u[vi0*3+2]);
		glVertex3d(v1[0]+u[vi1*3],v1[1]+u[vi1*3+1],v1[2]+u[vi1*3+2]);
	  }
	}
	glEnd();
	glEnable(GL_LIGHTING);
  }
}

void VolMeshDrawer::drawElement(const VolumetricMesh* volmesh){
  
  if(volmesh == NULL){
	WARN_LOG("the VolumetricMesh is null");
	return;
  }

  glDisable(GL_LIGHTING);
  glColor3d(0.8, 0.8, 0.8);
  glBegin(GL_LINES);
  const int num_ele = volmesh->numElements();
  const int num_ele_edges = volmesh->numElementEdges();
  int edge_index[num_ele_edges*2];
  for (int el = 0; el < num_ele; ++el){

	volmesh->getElementEdges(el,edge_index);
	for (int i = 0; i < num_ele_edges; ++i){
	  
	  const Vec3d &v0 = *(volmesh->vertex(edge_index[i*2]));
	  const Vec3d &v1 = *(volmesh->vertex(edge_index[i*2+1]));
	  glVertex3d(v0[0],v0[1],v0[2]);
	  glVertex3d(v1[0],v1[1],v1[2]);
	}
  }
  glEnd();
  glEnable(GL_LIGHTING);
}
