#ifndef _DRAGGROUPSELDRAW_H_
#define _DRAGGROUPSELDRAW_H_

#include <assertext.h>
#include <boost/foreach.hpp>
#include <volumetricMesh.h>
#include <vector>
#include <set>
using namespace std;

namespace LSW_SIM_UI{

  /**
   * @class DragGroupSelDraw privide "drawWithNames()" for a set of groups
   * constraint nodes for selecting constraint groups.
   * 
   */
  class DragGroupSelDraw{
	
  public:
	static void drawAllGroupsWithPoints(pVolumetricMesh_const vol_mesh,
										const vector<set<int> >&drag_groups, 
										const double*vol_u);

	static void drawAllGroupsWithShpere(pVolumetricMesh_const vol_mesh,
										const vector<set<int> >&drag_groups, 
										const double *vol_u, 
										const double radius);

  protected:
	static void drawPoints(pVolumetricMesh_const vol_mesh,
						   const set<int> &group,
						   const double *vol_u);

	static void drawSphere(pVolumetricMesh_const vol_mesh,
						   const set<int> &group,
						   const double *vol_u, const double radius);

	static void DrawSphere(float x, float y, float z, float fRadius, int M, int N);

  };
  
}//end of namespace

#endif /*_DRAGGROUPSELDRAW_H_*/
