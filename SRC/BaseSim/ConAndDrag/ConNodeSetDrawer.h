#ifndef _CONNODESETDRAWER_H_
#define _CONNODESETDRAWER_H_

#include <ConNodesOfFrameSet.h>

namespace LSW_SIM{
  
  /**
   * @class ConNodeSetDrawer render the the data of ConNodesOfFrameSet
   * 
   */
  class ConNodeSetDrawer{
	
  public:
	static void draw(const set<pConNodesOfFrame> &node_set);

	static void draw(const set<pConNodesOfFrame> &node_set, const float radius);
	
  };
    
}//end of namespace

#endif /*_CONNODESETDRAWER_H_*/
