#include "ObjBBox.h"
using namespace LSW_RENDERMESH;

BBox<float,Eigen::Vector3f> ObjBBox::getBBbox(const ObjRenderMesh &objmesh){
  
  BBox<float,Eigen::Vector3f> box;
  const float *vertex = objmesh.Vertices();
  const int v_num = objmesh.getVerticesNum();
  if(v_num > 0 && vertex != NULL){
  	box.reset(&vertex[3],v_num);
  }
  return box;
}

void ObjBBox::getBBbox(const ObjRenderMesh &objmesh,float min[3],float max[3]){

  const BBox<float,Eigen::Vector3f> box = getBBbox(objmesh);
  box.getMinConner(min);
  box.getMaxConner(max);
}
