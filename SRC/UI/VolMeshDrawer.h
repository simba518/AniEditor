#ifndef _VOLMESHDRAWER_H_
#define _VOLMESHDRAWER_H_

#include <vector>
#include <TetMesh.h>
using namespace std;
using namespace UTILITY;

namespace LSW_BASE_SIM{
  
  enum VOLMESHRENDERTYPE{POINT,EDGE,ELEMENT};

  /**
   * render TetMesh.
   */
  class VolMeshDrawer{
	
  public:
	static void render(pTetMesh_const volmesh,
					   VOLMESHRENDERTYPE t = POINT,
					   const double *u = NULL);

  protected:
	VolMeshDrawer(){}

	static void drawPoint(pTetMesh_const volmesh,const double *u);
	static void drawPoint(pTetMesh_const volmesh);

	static void drawEdge(pTetMesh_const volmesh,const double *u);
	static void drawEdge(pTetMesh_const volmesh);

	static void drawElement(pTetMesh_const volmesh,const double *u);
	static void drawElement(pTetMesh_const volmesh);

  };
  
}//end of namespace

#endif /*_VOLMESHDRAWER_H_*/
