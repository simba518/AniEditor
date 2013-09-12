#ifndef _VOLBBOX_H_
#define _VOLBBOX_H_

#include <boost/shared_ptr.hpp>
#include <BBox.h>
#include <eigen3/Eigen/Dense>
#include <volumetricMesh.h>
using UTILITY::BBox;

namespace LSW_RENDERMESH{
  
  class VolBBox{
	
  public:
	static boost::shared_ptr<VolBBox> getInstance(){
	  
	  if(p_instance == NULL){
		p_instance = boost::shared_ptr<VolBBox>(new VolBBox());
	  }
	  return p_instance;
	}

	BBox<double,Eigen::Vector3d > getBBbox(pVolumetricMesh_const vol_mesh)const;
	void getBBbox(pVolumetricMesh_const vol_mesh,double min[3],double max[3])const;
	
  protected:
    VolBBox(){}
	
  private:
	static boost::shared_ptr<VolBBox> p_instance;
  };
  
  typedef boost::shared_ptr<VolBBox> pVolBBox;
  
}//end of namespace

#endif /*_VOLBBOX_H_*/
