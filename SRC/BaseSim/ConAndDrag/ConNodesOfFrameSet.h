#ifndef _CONNODESOFFRAMESET_H_
#define _CONNODESOFFRAMESET_H_

#include <boost/shared_ptr.hpp>
#include <ConNodesOfFrame.h>

namespace LSW_SIM{
  
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

#endif /*_CONNODESOFFRAMESET_H_*/
