#ifndef _OBJBBOX_H_
#define _OBJBBOX_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include "ObjRenderMesh.h"
#include <BBox.h>
using UTILITY::BBox;

namespace LSW_RENDERMESH{

  /**
   * @class ObjBBox calculate the bounding box of an obj mesh.
   * 
   */
  class ObjBBox{
	
  public:
	static BBox<float,Eigen::Vector3f> getBBbox(const ObjRenderMesh &obj_mesh);
	static void getBBbox(const ObjRenderMesh &obj_mesh, float min[3], float max[3]);
	
  protected:
    ObjBBox(){}
	
  };
  
  typedef boost::shared_ptr<ObjBBox> pObjBBox;
  
}//end of namespace

#endif /*_OBJBBOX_H_*/
