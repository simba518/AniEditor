#include <iostream>
using namespace std;

#include <boost/foreach.hpp>
#include <assertext.h>
#include "BarycenterDragger.h"
using namespace LSW_SIM;

BarycenterDragger::BarycenterDragger(){
  
  mouse_start[0] = 0;
  mouse_start[1] = 0;
  mouse_start[2] = 0;
  rmDragGroup();
}

void BarycenterDragger::setConGroups(const vector<set<int> > &groups, const double *vol_u){

  this->groups = groups;
  u_c0 = computeBaryCenter(vol_u,groups);
  u_c = u_c0;
  rmDragGroup();
}

void BarycenterDragger::setDragGroup(const int group_id){

  assert_in(group_id, 0, (int)groups.size() - 1);
  sel_group_id = group_id;
}

/// after calling, the function hasDragedGroup() will return false.
void BarycenterDragger::rmDragGroup(){

  sel_group_id = -1;
}

void BarycenterDragger::startDrag(double m_x,double m_y,double m_z){

  mouse_start[0] = m_x;
  mouse_start[1] = m_y;
  mouse_start[2] = m_z;
  u_c = u_c0;
}

void BarycenterDragger::dragTo(double mouse_x,double mouse_y,double mouse_z){

  if(this->hasDragedGroup()){

	const Vector3d drag_distance=Vector3d(mouse_x,mouse_y,mouse_z)-mouse_start;
	u_c[sel_group_id] = u_c0[sel_group_id] + drag_distance;
  }
}

void BarycenterDragger::stopDrag(){

  u_c0 = u_c;
  rmDragGroup();
}

vector<Vector3d > BarycenterDragger::computeBaryCenter
(const double *vol_u, const vector<set<int> > &groups)const{
  
  vector<Vector3d > v;
  BOOST_FOREACH(const set<int> &g, groups){
	v.push_back( computeBaryCenter(vol_u,g) );
  }
  return v;
}

Vector3d BarycenterDragger::computeBaryCenter
(const double *vol_u, const set<int> &group)const{

  assert (vol_u != NULL);
  Vector3d barycenter(0,0,0);
  BOOST_FOREACH(int ele, group){
	barycenter[0] += vol_u[ele*3+0];
	barycenter[1] += vol_u[ele*3+1];
	barycenter[2] += vol_u[ele*3+2];
  }  

  const int ele_num = (int)group.size();
  assert (ele_num > 0);
  barycenter *= (1.0f/ele_num);
  
  return barycenter;
}

bool BarycenterDragger::hasDragedGroup()const{

  return (sel_group_id >=0 && sel_group_id < (int)u_c.size());
}

const set<int> &BarycenterDragger::getSelGroupNods()const{
  
  assert (hasDragedGroup());
  return groups[sel_group_id];
}

const Vector3d &BarycenterDragger::getTargetPos()const{

  if (hasDragedGroup()){
	return u_c[sel_group_id];
  }else{
	return mouse_start;
  }
}
