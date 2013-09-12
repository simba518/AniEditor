#ifndef _VOL2OBJ_H_
#define _VOL2OBJ_H_

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include <boost/shared_array.hpp>
#include <volumetricMesh.h>
#include <ObjMesh/ObjRenderMesh.h>

using namespace LSW_RENDERMESH;
using namespace std;

namespace LSW_SIM{
  
  ///interpolation from volume mesh to obj mesh
  class Vol2Obj{
	
  public:
	Vol2Obj(){
	  p_ObjRenderMesh = NULL;
	}
	void setObjMesh(const ObjRenderMesh *p_obj){
	  this->p_ObjRenderMesh = p_obj;
	}
	Vol2Obj(const ObjRenderMesh &m_ObjRenderMesh){
	  p_ObjRenderMesh = &m_ObjRenderMesh;
	}
	void setVolumeMesh(boost::shared_ptr< VolumetricMesh > _p_VolumetricMesh ){

	  assert(_p_VolumetricMesh != NULL);
	  p_VolumetricMesh = _p_VolumetricMesh;
	}
	bool precomputeObjMeshInterp();
	bool saveObjInterpDatas(string filename)const;
	bool loadObjInterpDatas(string filename);
	bool interpObjMesh(const double *u);
	const vector<float> &getuTarget()const{

	  return uTarget;
	}

  private:
	const ObjRenderMesh *p_ObjRenderMesh;
	boost::shared_ptr<VolumetricMesh > p_VolumetricMesh;
	boost::shared_array<int > interp_vertices;
	boost::shared_array<double > interp_weights;
	//interpolated displacements of obj vertices
	vector<float> uTarget; 
  };
  
  typedef boost::shared_ptr< Vol2Obj > pVol2Obj; 
  typedef boost::shared_ptr< const Vol2Obj > pVol2Obj_const; 

}//end of namespace

#endif /*_VOL2OBJ_H_*/
