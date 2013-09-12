#include <boost/foreach.hpp>
#include "BarycenterDraggerVec.h"
using namespace LSW_SIM;

void BarycenterDraggerVec::setConGroups
(const vector<set<int> > &groups, const VectorXd &vol_u){
  
  assert(checkValidDation(groups,vol_u));
  assert (vol_u.size() > 0);
  BarycenterDragger::setConGroups(groups,&vol_u[0]);
}

void BarycenterDraggerVec::stopDrag(){
  
  BarycenterDragger::stopDrag();
}

bool BarycenterDraggerVec::checkValidDation
(const vector<set<int> > &groups, const VectorXd &vol_u)const{

  BOOST_FOREACH(set<int> g, groups){
	BOOST_FOREACH(int ele, g){
	  if ( ele < 0 || ele*3+2 > vol_u.size() ){
		return false;
	  }
	}
  }
  return true;
}

VectorXd BarycenterDraggerVec::getUcVec()const{

  VectorXd uc_vec(u_c.size()*3);
  for (int i = 0; i < (int)u_c.size(); ++i){
	uc_vec[i*3+0] = u_c[i][0];
	uc_vec[i*3+1] = u_c[i][1];
	uc_vec[i*3+2] = u_c[i][2];
  }
  return uc_vec;
}
