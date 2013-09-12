#include <boost/filesystem.hpp>
#include "Vol2ObjAdaptor.h"
using namespace LSW_SIM;

bool Vol2ObjAdaptor::precompute(){

  bool succ = false;
  if(p_VolumetricMesh != NULL && p_ObjRenderMesh != NULL &&
	 p_ObjRenderMesh->getVerticesNum() > 0 ){

	int nTarget = p_ObjRenderMesh->getVerticesNum();
	const float * targetLocations = p_ObjRenderMesh->Vertices();
	int * temp_interp_vertices = NULL;
	double * temp_interp_weights = NULL;

	int numExternalVertices = 
	  p_VolumetricMesh->buildInterpolationWeights
	  (nTarget, targetLocations,
	   &(temp_interp_vertices),
	   &(temp_interp_weights),
	   -1, NULL, 1);

	interp_vertices = boost::shared_array<int> (temp_interp_vertices);
	interp_weights = boost::shared_array<double> (temp_interp_weights);

	succ = (numExternalVertices >= 0 && interp_weights != NULL&& interp_vertices != NULL);
  }
  return succ;
}

bool Vol2ObjAdaptor::saveObjInterpData(string filename)const{

  bool succ = false;
  if (checkValidation()){
	
	int nTarget = p_ObjRenderMesh->getVerticesNum();
	int numElementVertices_ = p_VolumetricMesh->numElementVertices();
  
	const int return_code = VolumetricMesh::saveInterpolationWeights
	  (filename.c_str(), nTarget, numElementVertices_, 
	   interp_vertices.get(), interp_weights.get());
	succ = (return_code == 0);
  }
  return succ;
}

bool Vol2ObjAdaptor::loadObjInterpData(string filename){

  bool succ = false;
  if( boost::filesystem::exists(filename) && checkValidation()){

	const int nTarget = p_ObjRenderMesh->getVerticesNum();
	const int numElementVertices_ = p_VolumetricMesh->numElementVertices();

	int * temp_interp_vertices = NULL;
	double * temp_interp_weights = NULL;

	const int code = VolumetricMesh::loadInterpolationWeights
	  (filename.c_str(), nTarget, numElementVertices_, 
	   &temp_interp_vertices, &temp_interp_weights);

	interp_vertices = boost::shared_array<int> (temp_interp_vertices);
	interp_weights = boost::shared_array<double> (temp_interp_weights);

	succ =  (code == 0);	
  }else{
	printf("error:in loadObjInterpDatas(char *),file %s is not exist\n",filename.c_str());
  }
  return succ;
}

bool Vol2ObjAdaptor::interpolate(const double *vol_u,pObjRenderMesh obj_mesh){

  vector<float> uTarget;
  bool succ = false;
  if(interpolate(vol_u,uTarget) && obj_mesh != NULL){
  	succ = obj_mesh->move(uTarget);
  }
  return succ;
}

bool Vol2ObjAdaptor::interpolate(const double *vol_u,vector<float> &uTarget){

  if(p_ObjRenderMesh != NULL){
	const int nTarget = p_ObjRenderMesh->getVerticesNum();
	uTarget.resize(nTarget*3);	
  }
  return interpolate(vol_u,&(uTarget[0]));
}

bool Vol2ObjAdaptor::interpolate(const double *vol_u,float* uTarget){

  bool succ = false;
  //interp_weights and interp_vertices should be computed already.
  if(checkValidation() && 
	 interp_weights != NULL && interp_vertices != NULL && 
	 uTarget != NULL && vol_u!=NULL){

	const int nTarget = p_ObjRenderMesh->getVerticesNum();	
	const int numElementVertices_ = p_VolumetricMesh->numElementVertices();
	VolumetricMesh::interpolate(vol_u, &uTarget[0], nTarget, 
								numElementVertices_, interp_vertices.get(), 
								interp_weights.get());
	succ = true;
  }
  return succ;
}

bool Vol2ObjAdaptor::checkValidation()const{

  return (p_VolumetricMesh != NULL)&&(p_ObjRenderMesh != NULL)&&
	(p_ObjRenderMesh->getVerticesNum() > 0);
}
