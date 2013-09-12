#ifndef _VOLMESHDRAWER_H_
#define _VOLMESHDRAWER_H_

#include <vector>
#include <volumetricMesh.h>
using namespace std;

namespace LSW_BASE_SIM{
  
  enum VOLMESHRENDERTYPE{POINT,EDGE,ELEMENT};

  /**
   * render VolumetricMesh.
   */
  class VolMeshDrawer{
	
  public:
	static void render(const VolumetricMesh* volmesh,
					   VOLMESHRENDERTYPE t = POINT,
					   const double *u = NULL);

  protected:
	VolMeshDrawer(){}

	static void drawPoint(const VolumetricMesh* volmesh,const double *u);
	static void drawPoint(const VolumetricMesh* volmesh);

	static void drawEdge(const VolumetricMesh* volmesh,const double *u);
	static void drawEdge(const VolumetricMesh* volmesh);

	static void drawElement(const VolumetricMesh* volmesh,const double *u);
	static void drawElement(const VolumetricMesh* volmesh);

  };
  
}//end of namespace

#endif /*_VOLMESHDRAWER_H_*/
