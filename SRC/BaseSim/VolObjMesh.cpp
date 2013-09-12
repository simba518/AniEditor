#include <JsonFilePaser.h>
#include "VolObjMesh.h"
using namespace LSW_SIM;

bool VolObjMesh::initialize(const string init_filename){
  
  UTILITY::JsonFilePaser inf;
  bool succ = false;
  if ( inf.open(init_filename) ){

	string filename;
	if ( inf.readFilePath("obj_mesh_file",filename,true) ){
	  succ = loadObjMesh(filename);
	}
	if ( inf.readFilePath("vol_filename",filename,true) ){
	  succ &= loadVolMesh(filename);
	}
	if ( inf.readFilePath("vol_surface",filename,true) ){
	  succ &= loadVolSurface(filename);
	}
	if ( inf.readFilePath("vol2obj_weights_filename",filename,true)){
	  succ &= loadInterpWeights(filename);
	}
  }
  return succ;
}

int VolObjMesh::objVertexNum()const{
 
  if (obj_mesh){
	return obj_mesh->getVerticesNum();
  }
  return 0;
}

int VolObjMesh::volVertexNum()const{

  if(vol_mesh){
	return vol_mesh->numVertices();
  }
  return 0;
}

int VolObjMesh::volElementNum()const{
  
  if(vol_mesh){
	return vol_mesh->numElements();
  }
  return 0;
}

double VolObjMesh::getMaxRadius()const{

  double r = max[0] - min[0];
  r = r >= (max[1] - min[1]) ? r : (max[1] - min[1]);
  r = r >= (max[2] - min[2]) ? r : (max[2] - min[2]);
  return r >= 0 ? r : 0;
}
