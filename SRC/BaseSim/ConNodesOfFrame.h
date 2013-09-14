#ifndef _CONNODESOFFRAME_H_
#define _CONNODESOFFRAME_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <SelectionGroup.h>
using namespace Eigen;
using namespace UTILITY;

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


  /**
   * @class ConNodesOfFrameSet set of ConNodesOfFrame
   * 
   */
  class ConNodesOfFrameSet{
	
  public:
	ConNodesOfFrameSet(){
	  
	}
	ConNodesOfFrameSet(const ConNodesOfFrameSet &other);
	ConNodesOfFrameSet &operator=(const ConNodesOfFrameSet &other);
	ConNodesOfFrameSet &copy(const ConNodesOfFrameSet &other);
	void clear(){
	  node_set.clear();
	}
	
	bool load(const string filename);
	bool save(const string filename)const;
	bool saveConPositions(const string filename)const;

	void addConNodeGroup(pConNodesOfFrame node_group);
	void addConNodeGroup(const ConNodesOfFrame &node_group);
	void addConNodeGroup( const vector<int> &con_nodes, const int frame_id );
	void rmConNodeGroup( const vector<int> &con_nodes, const int frame_id );
	void setConPos(const VectorXd &uc, const VectorXd &rest,const int frame_id);
	bool updateConNodeGroup(const VectorXd &uc,const int frame_id);

	pConNodesOfFrame_const getConNodeGroup(const int frame_id)const;
	pConNodesOfFrame getConNodeGroup(const int frame_id);
	const set<pConNodesOfFrame> &getConNodeGroups()const{
	  return node_set;
	}
	vector<ConNodesOfFrame> getSortedConNodeGroups()const;
	
	bool equal (const ConNodesOfFrameSet &other,const double tol = 1e-6)const;

  private:
	set<pConNodesOfFrame> node_set;
  };
  
}//end of namespace

#endif /*_CONNODESOFFRAME_H_*/
