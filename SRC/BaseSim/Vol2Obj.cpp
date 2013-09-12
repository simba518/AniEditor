#include <boost/filesystem.hpp>
#include <Log.h>
#include "Vol2Obj.h"
using namespace LSW_SIM;

bool Vol2Obj::precomputeObjMeshInterp(){
  
  bool succ = false;
  if(p_VolumetricMesh != NULL && p_ObjRenderMesh != NULL &&
	 p_ObjRenderMesh->getVerticesNum() > 0 ){

	int nTarget = p_ObjRenderMesh->getVerticesNum();
	const float * targetLocations = p_ObjRenderMesh->Vertices();
	int * temp_interp_vertices = NULL;
	double * temp_interp_weights = NULL;

	const int numExternalVertices = p_VolumetricMesh->buildInterpolationWeights
	  (nTarget, targetLocations,&(temp_interp_vertices),&(temp_interp_weights),-1, NULL, 1);

	interp_vertices = boost::shared_array<int> (temp_interp_vertices);
	interp_weights = boost::shared_array<double> (temp_interp_weights);

	succ = (numExternalVertices >= 0 && interp_weights != NULL && interp_vertices != NULL);
  }
  return succ;
}

///save interp_weights interp_vertices
bool Vol2Obj::saveObjInterpDatas(string filename)const{

  assert(p_VolumetricMesh != NULL);
  assert(p_ObjRenderMesh != NULL);
  assert_gt(p_ObjRenderMesh->getVerticesNum(), 0);
   
  int nTarget = p_ObjRenderMesh->getVerticesNum();
  int numElementVertices_ = p_VolumetricMesh->numElementVertices();
  
  const int return_code = VolumetricMesh::saveInterpolationWeights
	(filename.c_str(), nTarget, numElementVertices_, 
	 interp_vertices.get(), interp_weights.get());

  return (return_code == 0);
}

///read interpolation weights from the disk
bool Vol2Obj::loadObjInterpDatas(string filename){
  
  TRACE_FUN();

  assert(p_VolumetricMesh != NULL);
  assert(p_ObjRenderMesh != NULL);

  bool succ = false;
  if(boost::filesystem::exists(filename)){

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
	ERROR_LOG("file"<<filename <<" not found\n");
  }
  return succ;
}

/** 
 * interpolate displacement from VolumetricMesh vertices to obj mesh.
 * @param[in] u displacement of VolumetricMesh vertices.
 * @return the vertices of m_ObjRenderMesh will be updated.
 */
bool Vol2Obj::interpObjMesh(const double *u){

  bool succ = false;
  if(p_VolumetricMesh != NULL && p_ObjRenderMesh != NULL){

	//interp_weights and interp_vertices should be computed already.
	if( !(p_ObjRenderMesh->getVerticesNum() <=0 || interp_weights == NULL || interp_vertices == NULL) ){
	  
	  const int numElementVertices_ = p_VolumetricMesh->numElementVertices();
	  const int nTarget = p_ObjRenderMesh->getVerticesNum();
	  uTarget.resize(nTarget*3);
	  VolumetricMesh::interpolate(u, &uTarget[0], nTarget, 
								  numElementVertices_, interp_vertices.get(), 
								  interp_weights.get());
	  succ = true;
	}
  }
  return succ;
}
