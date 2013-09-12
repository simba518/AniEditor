#include <iostream>
#include <memory.h>
#include <PhongShader.h>
#include <Log.h>
#include <assertext.h>
#include "ObjRenderMesh.h"
using namespace std;
using namespace LSW_RENDERMESH;

ObjRenderMesh::ObjRenderMesh(){

  obj_model = NULL;
  default_render_mode = GLM_NONE;
  use_phone_shading = true;
}

ObjRenderMesh::ObjRenderMesh(const string &file_name){

  read(file_name);
}

bool ObjRenderMesh::read(const string &file_name){

  TRACE_FUN();
  this->file_path = file_name;

  clean();
  obj_model = glmReadOBJ(const_cast<char* >(file_name.c_str()));
  bool succ = true;
  if (obj_model->numvertices <= 0){
	clean();
	succ = false;
  }
  if ( succ ){
  	//compute normals
  	if (obj_model->numfacetnorms <= 0){
  	  glmFacetNormals(obj_model);
  	}
  	if (obj_model->numnormals <= 0){
  	  glmVertexNormals(obj_model,90,false);
  	}
	// init original positions of the vertices
	ori_pos.resize(obj_model->numvertices*3);
	for(size_t i = 0;i < obj_model->numvertices; i++){
	  ori_pos[i*3+0] = obj_model->vertices[i*3+3+0];
	  ori_pos[i*3+1] = obj_model->vertices[i*3+3+1];
	  ori_pos[i*3+2] = obj_model->vertices[i*3+3+2];
	}
	// init defualt render mode
	initDefaultRenderMode();
  }
  return succ;
}

bool ObjRenderMesh::save(const string &file_name,int mode)const{
  
  if(obj_model != NULL){
	glmWriteOBJ(obj_model,const_cast<char*>(file_name.c_str()),mode);
	return true;
  }
  return false;
}

float ObjRenderMesh::unitize(){

  if(obj_model != NULL){
	return glmUnitize(obj_model);
  }
  return 1.0f;
}

void ObjRenderMesh::scale(const float scale_f){

  if(obj_model != NULL){
	glmScale(obj_model,scale_f);
  }
}

void ObjRenderMesh::getDimension(float dimensions[3])const{

  if(obj_model != NULL){
	glmDimensions(obj_model,dimensions);
  }else{
	dimensions[0] = 0;
	dimensions[1] = 0;
	dimensions[2] = 0;
  }
}

void ObjRenderMesh::computeNormals(){

  if(obj_model != NULL){
	glmFacetNormals(obj_model);
	glmVertexNormals(obj_model,90.0f,false);	
  }
}

/** 
 * reset the obj model's vertices: vertices = ori_pos + u.
 * @param u the displacements of the vertices.
 * @return the parameter (dimension) is valid or not.
 *
 * @note ori_pos is initialized at initialize() or ObjRenderMesh<float>::read()
 */
bool ObjRenderMesh::move(const vector<float> &u){

  assert_eq((int)u.size(), this->getVerticesNum()*3 );
  assert_eq(ori_pos.size(), u.size());
	  
  for ( size_t i = 0; i < u.size(); i++ ){
	// yes, vertices start from [3]
  	obj_model->vertices[i+3] = ori_pos[i]+u[i];
  }
  return true;
}

bool ObjRenderMesh::setVertices(const float *v,int len){
  
  assert(v != NULL);
  assert_eq(len, this->getVerticesNum()*3 );

  for ( int i = 0; i < len; i++ ){
	// yes, vertices start from [3]
  	obj_model->vertices[i+3] = v[i];
  }
  return true;
}

bool ObjRenderMesh::setVertices(const vector<float> &vertices){

  assert_gt(vertices.size(),0);
  return setVertices(&vertices[0],vertices.size());
}

void ObjRenderMesh::reset(){

  if(obj_model != NULL){

	assert_eq((int)ori_pos.size(), this->getVerticesNum()*3);
	for ( size_t i = 0; i < ori_pos.size(); i++ ){
	  obj_model->vertices[i+3] = ori_pos[i];
	}	
  }
}

const float *ObjRenderMesh::Vertices()const{

  if(obj_model == NULL || obj_model->numvertices <=0){
  	return NULL;
  }else{
	//yes, it start from numvertices[3]
	return &(obj_model->vertices[3]);
  }
}

void ObjRenderMesh::clean(){

  if(obj_model != NULL){
	glmDelete(obj_model);
  }
  obj_model = NULL;
  ori_pos.clear();
  default_render_mode = GLM_NONE;
}

void ObjRenderMesh::draw(const int mode)const{

  if(getVerticesNum() > 0){

	if (usePhoneShading()){

	  PhongShader::getInstance()->bind();
	  glmDraw(obj_model,mode);
	  PhongShader::getInstance()->unbind();

	}else{
	  glmDraw(obj_model,mode);
	}
  }
}

void ObjRenderMesh::draw()const{
  
  draw(default_render_mode);
}

void ObjRenderMesh::draw(float color_r,float color_g,float color_b)const{

  assert(obj_model != NULL);
  GLMmaterial *mtl = obj_model->materials;
  const int num_mtl = obj_model->nummaterials;
  assert_ge(num_mtl,0);
  vector<float> diffuse(num_mtl*4);
  vector<float> ambient(num_mtl*4);
  for (int i = 0; i < num_mtl; ++i){

	diffuse[i*4+0] = mtl[i].diffuse[0];
	diffuse[i*4+1] = mtl[i].diffuse[1];
	diffuse[i*4+2] = mtl[i].diffuse[2];
	diffuse[i*4+3] = mtl[i].diffuse[3];

	ambient[i*4+0] = mtl[i].ambient[0];
	ambient[i*4+1] = mtl[i].ambient[1];
	ambient[i*4+2] = mtl[i].ambient[2];
	ambient[i*4+3] = mtl[i].ambient[3];

	mtl[i].diffuse[0] = color_r;
	mtl[i].diffuse[1] = color_g;
	mtl[i].diffuse[2] = color_b;
	mtl[i].diffuse[3] = 1.0f;

	mtl[i].ambient[0] = color_r;
	mtl[i].ambient[1] = color_g;
	mtl[i].ambient[2] = color_b;
	mtl[i].ambient[3] = 1.0f;
  }

  glmDraw(obj_model,default_render_mode|GLM_MATERIAL);

  for (int i = 0; i < num_mtl; ++i){

	mtl[i].diffuse[0] = diffuse[i*4+0];
	mtl[i].diffuse[1] = diffuse[i*4+1];
	mtl[i].diffuse[2] = diffuse[i*4+2];
	mtl[i].diffuse[3] = diffuse[i*4+3];

	mtl[i].ambient[0] = ambient[i*4+0];
	mtl[i].ambient[1] = ambient[i*4+1];
	mtl[i].ambient[2] = ambient[i*4+2];
	mtl[i].ambient[3] = ambient[i*4+3];
  }
}

int ObjRenderMesh::getVerticesNum()const{
  
  if(obj_model != NULL){
	return obj_model->numvertices;
  }
  return 0;
} 

int ObjRenderMesh::getTriangleNum()const{
  
  if(obj_model != NULL){
	return obj_model->numtriangles;
  }
  return 0;  
} 

int ObjRenderMesh::getMaterialNum()const{
  
  if(obj_model != NULL){
	return obj_model->nummaterials;
  }
  return 0;
} 

int ObjRenderMesh::getTextureNum()const{
  
  if(obj_model != NULL){
	return obj_model->numtexcoords;
  }
  return 0;
}  

void ObjRenderMesh::initDefaultRenderMode(){
  
  default_render_mode = GLM_NONE;
  if(getVerticesNum() > 0){
	default_render_mode = (default_render_mode | GLM_SMOOTH );
  }
  if(getMaterialNum() > 0){
	default_render_mode = (default_render_mode | GLM_MATERIAL);
  }else{
	default_render_mode = (default_render_mode | GLM_COLOR);
  }
  if(getTextureNum() > 0){
	default_render_mode = (default_render_mode | GLM_TEXTURE);
  }
}

ObjRenderMesh::~ObjRenderMesh(){
  
  this->clean();
}

void ObjRenderMesh::getCenter(double center_pos[3])const{
  
  const float *v = this->Vertices();
  center_pos[0] = 0;
  center_pos[1] = 0;
  center_pos[2] = 0;
  if (v != NULL){

	const int n = getVerticesNum();
	for (int i = 0; i < n; ++i){

	  center_pos[0] += v[i*3+0];
	  center_pos[1] += v[i*3+1];
	  center_pos[2] += v[i*3+2];
	}	

	if (n > 0){

	  center_pos[0] /= double(n);
	  center_pos[1] /= double(n);
	  center_pos[2] /= double(n);
	}
  }
}

void ObjRenderMesh::moveCenterTo(const double center_pos[3]){
  
  scale( unitize() ); // move to original

  // move to center_pos
  float *v = &(obj_model->vertices[3]);
  if (v != NULL){

	const int n = getVerticesNum();
	for (int i = 0; i < n; ++i){
	  v[i*3+0] += center_pos[0];
	  v[i*3+1] += center_pos[1];
	  v[i*3+2] += center_pos[2];
	}
  }
}

void ObjRenderMesh::Triangles(vector<int> &tri)const{
  
  tri.resize(0);
  if(obj_model){
	tri.resize(getTriangleNum()*3);
	vector<GLuint> faces(getTriangleNum()*3);
	// @todo glmGetFaces is defined by lsw in glm.
	glmGetFaces((GLMmodel*)obj_model,&faces[0]);
	for (size_t i = 0; i < faces.size(); ++i){
	  tri[i] = (int)faces[i]-1;
	}
  }
}
