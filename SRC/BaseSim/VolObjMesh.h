#ifndef _VOLOBJMESH_H_
#define _VOLOBJMESH_H_

#include <boost/shared_ptr.hpp>
#include <volumetricMeshLoader.h>
#include <VolBBox.h>
#include <Vol2Obj.h>

namespace LSW_SIM{
  
  /**
   * @class VolObjMesh contain volume mesh ,surface mesh and obj mesh, and the
   * interpolate mehtod. The number and index of the surface mesh is equal to
   * the volumetric mesh, and can be obtained by using the Barbic's application.
   * 
   */
  class VolObjMesh{
	
  public: 
	bool initialize(const string filename);
	bool loadObjMesh(const string filename){
	  if (obj_mesh == NULL)
		obj_mesh = pObjRenderMesh(new ObjRenderMesh() );
	  return obj_mesh->read(filename);
	}
	bool loadVolSurface(const string filename){
	  if (vol_surface == NULL)
		vol_surface = pObjRenderMesh(new ObjRenderMesh() );
	  return vol_surface->read(filename);
	}
	bool loadVolMesh(const string filename){
	  vol_mesh = pVolumetricMesh(VolumetricMeshLoader::load(filename.c_str()));
	  if (vol_mesh)
		VolBBox::getInstance()->getBBbox(vol_mesh,min,max);
	  return (vol_mesh != NULL);
	}
	bool loadInterpWeights(const string filename){
	  if(initVol2Obj())
		return vol2obj.loadObjInterpDatas(filename);
	  return false;
	}
	bool saveInterpWeights(const string filename)const{
	  return vol2obj.saveObjInterpDatas(filename);
	}

	pVolumetricMesh getVolMesh()const{
	  return vol_mesh;
	}
	pObjRenderMesh getObjMesh()const{
	  return obj_mesh;
	}
	pObjRenderMesh getVolSurface()const{
	  return vol_surface;
	}

	bool precomputeInterpWeights(){
	  if(initVol2Obj())
		return vol2obj.precomputeObjMeshInterp();
	  return false;
	}
	bool interpolate(const double *vol_u){
	  if ( vol2obj.interpObjMesh(vol_u) )
		return obj_mesh->move(vol2obj.getuTarget());
	  return false; 
	}

	int objVertexNum()const;
	int volVertexNum()const;
	int volElementNum()const;

	void getBBbox(double min[3],double max[3])const{
	  min[0] = this->min[0];  min[0] = this->min[1];  min[0] = this->min[2];  
	  max[0] = this->max[0];  max[0] = this->max[1];  max[0] = this->max[2];
	}
	double getMaxRadius()const;

  protected:
	bool initVol2Obj(){
	  
	  if(obj_mesh != NULL && vol_mesh != NULL){
		vol2obj.setObjMesh(obj_mesh.get());
		vol2obj.setVolumeMesh(vol_mesh);
		return true;
	  }		
	  return false;
	}
	
  private:
	pObjRenderMesh obj_mesh;
	pObjRenderMesh vol_surface;
	pVolumetricMesh vol_mesh;
	Vol2Obj vol2obj;

	// bound box, initialized when volmesh loaded.
	double min[3]; // min value for x,y,z
	double max[3]; // max value for x,y,z
  };
  
  typedef boost::shared_ptr<VolObjMesh> pVolObjMesh;
  typedef boost::shared_ptr<const VolObjMesh> pVolObjMesh_const;
  
}//end of namespace

#endif /*_VOLOBJMESH_H_*/
