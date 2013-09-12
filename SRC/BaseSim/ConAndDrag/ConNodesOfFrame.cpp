#include <assertext.h>
#include "ConNodesOfFrame.h"
using namespace LSW_SIM;

void ConNodesOfFrame::setConNodes
(const vector<set<int> > &con_nodes_set, const VectorXd &b_uc, 
 const VectorXd &rest_barycenters,const int frame_id){
  
  assert_ge (frame_id,0);
  assert_eq (b_uc.size(), (int)(con_nodes_set.size()*3));
  assert_eq (rest_barycenters.size(), b_uc.size());
  
  this->con_nodes_set.setGroup( con_nodes_set );
  this->barycenter_uc = b_uc;
  this->barycenter_rest = rest_barycenters;
  this->frame_id = frame_id;
}

VectorXd ConNodesOfFrame::getBarycenter()const{

  return barycenter_rest + barycenter_uc;
}

bool ConNodesOfFrame::equal(const ConNodesOfFrame &other,const double tol)const{

  const vector<set<int> > &other_con_nodes_set = other.getConNodesSet();
  const VectorXd &other_barycenter_uc = other.getBarycenterUc();
  const VectorXd &other_barycenter = other.getBarycenter();
  const VectorXd this_barycenter = this->getBarycenter();
  const int other_frame_id = other.getFrameId();

  const bool is_equal =  (other_frame_id == frame_id)&&
	((other_barycenter-this_barycenter).norm()<tol)&&
	((other_barycenter_uc-barycenter_uc).norm()<tol)&&
	(other_con_nodes_set == this->getConNodesSet());

  return is_equal;
}
