#ifndef _BARYCENTERDRAGGER_H_
#define _BARYCENTERDRAGGER_H_

#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;

namespace LSW_SIM{
  
  /**
   * @class BarycenterDragger drag the barycenter of one of the selected
   * constrained groups.
   * 
   */
  class BarycenterDragger{
	
  public:
	BarycenterDragger();
	void setConGroups(const vector<set<int> > &groups, const double *vol_u);
	void setDragGroup(const int group_id);
	void rmDragGroup();
	void startDrag(double mouse_x,double mouse_y,double mouse_z);
	void dragTo(double mouse_x,double mouse_y,double mouse_z);
	void stopDrag();
	
	const vector<Vector3d> &getUc()const{
	  return u_c;
	}
	int getSelGroup()const{
	  return sel_group_id;
	}
	const set<int> &getSelGroupNods()const;
	void removeAllConGroups(){
	  groups.clear();
	  sel_group_id = -1;
	}
	bool hasDragedGroup()const;
	const Vector3d &getTargetPos()const;
	
  protected:
	vector<Vector3d > computeBaryCenter(const double *vol_u, const vector<set<int> > &groups)const;
	Vector3d computeBaryCenter(const double *vol_u, const set<int> &group)const;
	
  protected:
	vector<Vector3d > u_c0;  // the barycenters of the displacements of the
							  // constrained groups before draging.

	vector<Vector3d > u_c;   // the barycenters of the displacements of the
							  // constrained groups after draging.

	Vector3d mouse_start;    // the mouse position in the world coordinates
							  // when drag begin.

	vector<set<int> > groups; // all of the constrained groups.
	int sel_group_id; // selected group for draging.
	
  };
  
  typedef boost::shared_ptr<BarycenterDragger> pBarycenterDragger;
  
}//end of namespace

#endif /*_BARYCENTERDRAGGER_H_*/
