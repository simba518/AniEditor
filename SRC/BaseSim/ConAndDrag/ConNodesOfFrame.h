#ifndef _CONNODESOFFRAME_H_
#define _CONNODESOFFRAME_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
using namespace Eigen;

#include <SelectionGroup.h>
using namespace LSW_SIM;

namespace LSW_SIM{
  
  /**
   * @class ConNodesOfFrame record the constrained nodes and target
   * displacements of the barycenter of one frame.
   * 
   */
  class ConNodesOfFrame{
	
  public:
	ConNodesOfFrame(const int frame_id = 0){
	  this->frame_id = frame_id;
	}
	ConNodesOfFrame
	(const vector<set<int> > &con_nodes_set, const VectorXd &barycenter_uc,
	 const VectorXd &rest_barycenters, const int frame_id){
	  setConNodes(con_nodes_set,barycenter_uc,rest_barycenters,frame_id);
	}

	void setConNodes
	(const vector<set<int> > &con_nodes_set,const VectorXd &barycenter_uc,
	 const VectorXd &rest_barycenters,const int f);

	void addConNodeGroup( const vector<int> &con_nodes){
	  con_nodes_set.addGroup(con_nodes);
	}
	void rmConNodeGroup( const vector<int> &con_nodes){
	  con_nodes_set.removeGroup(con_nodes);
	}
	void updateConPos(const VectorXd &bary_uc,const VectorXd &bary_rest){
	  setUc(bary_uc);
	  this->barycenter_rest = bary_rest;
	}
	void setUc(const VectorXd &barycenter_uc){
	  this->barycenter_uc = barycenter_uc;
	}

	void clear(){
	  con_nodes_set.clear();
	  barycenter_uc.resize(0);
	  frame_id = 0;
	}

	bool isEmpty()const{
	  return con_nodes_set.isEmpty();
	}
	const vector<set<int> > &getConNodesSet()const{
	  return con_nodes_set.getGroup();
	}
	const VectorXd &getBarycenterUc()const{
	  return barycenter_uc;
	}
	const VectorXd &getBarycenterRest()const{
	  return barycenter_rest;
	}
	VectorXd getBarycenter()const;
	const int getFrameId()const{
	  return frame_id;
	}

	bool equal (const ConNodesOfFrame &other,const double tol = 1e-6)const;
	
	bool operator< (const ConNodesOfFrame &other)const{
	  return this->getFrameId() < other.getFrameId();
	}

  private:
	// Each set of constrained nodes record one group of constrained nodes. The
	// i-th 3x1 subvector in barycenter_uc represents the barycenter of the
	// displacements of i-th group of con_nodes_set (to the rest shape).
	SelectionGroup con_nodes_set;
	VectorXd barycenter_uc;
	VectorXd barycenter_rest; // the barycenteres of connodes of the rest shape.
	int frame_id;   // frame index
  };
  
  typedef boost::shared_ptr<ConNodesOfFrame> pConNodesOfFrame;
  typedef boost::shared_ptr<const ConNodesOfFrame> pConNodesOfFrame_const;
  
}//end of namespace

#endif /*_CONNODESOFFRAME_H_*/
