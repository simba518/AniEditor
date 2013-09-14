#ifndef _VOL2OBJADAPTOR_H_
#define _VOL2OBJADAPTOR_H_

#include <vector>
#include <string>
#include <volumetricMesh.h>
#include <ObjRenderMesh.h>
#include <boost/shared_array.hpp>

using namespace LSW_RENDERMESH;
using namespace std;

namespace LSW_SIM{
  
  /**
   * @class Vol2ObjAdaptor a class interpolating a objmesh from a volume mesh,
   * similar to the Vol2Obj class, but with much elagent interfaces.
   */
  class Vol2ObjAdaptor{
	
  public:
	Vol2ObjAdaptor(){}
	Vol2ObjAdaptor(pObjRenderMesh_const obj_mesh,pVolumetricMesh vol_mesh){
	  setMesh(obj_mesh,vol_mesh);
	}
	void setMesh(pObjRenderMesh_const obj_mesh,pVolumetricMesh vol_mesh){
	  setObjMesh(obj_mesh);
	  setVolumeMesh(vol_mesh);
	}
	void setObjMesh(pObjRenderMesh_const obj_mesh){
	  this->p_ObjRenderMesh = obj_mesh;
	}
	void setVolumeMesh(pVolumetricMesh vol_mesh){
	  this->p_VolumetricMesh = vol_mesh;
	}
	bool precompute();
	bool saveObjInterpData(string filename)const;
	bool loadObjInterpData(string filename);
	bool interpolate(const double *vol_u,pObjRenderMesh obj_mesh);
	bool interpolate(const double *vol_u,vector<float> &obj_u);
	bool interpolate(const double *vol_u,float* obj_u);
	
  protected:
	bool checkValidation()const;

  private:
	pObjRenderMesh_const p_ObjRenderMesh;
	pVolumetricMesh p_VolumetricMesh;
	boost::shared_array<int > interp_vertices;
	boost::shared_array<double > interp_weights;
  };
  
  typedef boost::shared_ptr<Vol2ObjAdaptor> pVol2ObjAdaptor;
  typedef boost::shared_ptr<const Vol2ObjAdaptor> pVol2ObjAdaptor_const;
  
}//end of namespace

#endif /*_VOL2OBJADAPTOR_H_*/
