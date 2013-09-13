#include <Log.h>
#include "RSWarperExt.h"
#include "RSCoordComp.h"
using namespace IEDS;

void RSWarperExt::setConNodes(const vector<set<int> > &c_nodes,const VectorXd &bcen_uc){
  
  rs2euler.setConNodes(c_nodes, bcen_uc);
  rs2euler.precompute();
}

void RSWarperExt::setFixedNodes(const vector<int>&c_nodes){

  rs2euler.setFixedNodes(c_nodes);
}

void RSWarperExt::setRefU(const vector<VectorXd> &ref_u){
  
  this->ref_y = ref_u;
}
	
bool RSWarperExt::precompute(){
  
  TRACE_FUN();
  bool succ = rs2euler.precompute();
  if (succ){
	for (int i = 0; i < (int)ref_y.size(); ++i){
	  const VectorXd u = ref_y[i];
	  computeRSCoord( u , ref_y[i] , false);
	}
  }
  return succ;
}

void RSWarperExt::warp(const VectorXd &p,int frame_id, VectorXd &u){

  VectorXd y;
  LSW_WARPING::RSCoordComp::constructWithWarp(rs2euler.get_G(),p,y);

  assert_in (frame_id, 0, (int)ref_y.size()-1);
  assert_eq(y.size(),ref_y[frame_id].size());
  y += ref_y[frame_id];

  rs2euler.reconstruct(y,u);
  assert_eq(u.size(),p.size());
}

/** 
 * compute the RS coordinates of u.
 * 
 * @param u the displacement in full space with respect to the rest shape.
 * @param RS_u the RS-coordinates compute from u.
 */
void RSWarperExt::computeRSCoord(const VectorXd &u, VectorXd &RS_u, const bool linear_disp)const{
  
  const SparseMatrix<double> &G = rs2euler.get_G();
  if (linear_disp){
	LSW_WARPING::RSCoordComp::constructWithWarp(G,u,RS_u);
  }else{
	LSW_WARPING::RSCoordComp::constructWithoutWarp(G,u,RS_u);
  }
}

const VectorXd &RSWarperExt::getInputRSCoord(const int frame_id)const{
  
  assert_in (frame_id,0,(int)ref_y.size()-1);
  return ref_y[frame_id];
}
