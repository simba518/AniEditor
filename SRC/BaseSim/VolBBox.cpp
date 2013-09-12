#include <iostream>
#include "VolBBox.h"
using namespace std;
using namespace LSW_RENDERMESH;

pVolBBox VolBBox::p_instance;

BBox<double,Eigen::Vector3d> VolBBox::getBBbox(pVolumetricMesh_const vol_mesh)const{

  BBox<double,Eigen::Vector3d> box;
  if(vol_mesh != NULL && vol_mesh->numVertices() > 0){

	vector<double> vertices;
	vol_mesh->vertices(vertices);
	box.reset(vertices);
  }
  return box;  
}

void VolBBox::getBBbox(pVolumetricMesh_const vol_mesh,double min[3],double max[3])const{

  const BBox<double,Eigen::Vector3d> box = getBBbox(vol_mesh);
  box.getMinConner(min);
  box.getMaxConner(max);  
}
