#ifndef _DRAGTRAJECTORYRECORD_H_
#define _DRAGTRAJECTORYRECORD_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include <assertext.h>
using namespace std;
using namespace Eigen;

namespace LSW_ANI_EDITOR{

  class OneDragRecord{
	
  public:
    OneDragRecord(){
	  frame_num = 0;
	}
	OneDragRecord(const int frame_num, const VectorXd &uc){
	  this->frame_num = frame_num;
	  this->uc = uc;
	}
	int getFrameNum()const{
	  return frame_num;
	}
	const VectorXd &getUc()const{
	  return uc;
	}
	
  private:
	int frame_num;
	VectorXd uc;
  };

  /**
   * @class DragTrajectoryRecord Records the dragging trajectory as pairs as:
   * frame-number and uc, and provides the interfaces for saving and loading the
   * records.
   * 
   */
  class DragTrajectoryRecord{
	
  public:
	DragTrajectoryRecord(){}
	void record(const int frame_num, const VectorXd &uc){
	  records.push_back(OneDragRecord(frame_num, uc));
	}
	const OneDragRecord &getRecord(const int record_num)const{
	  assert_in (record_num,0,totalRecord()-1);
	  return records[record_num];
	}
	int totalRecord()const{
	  return (int)records.size();
	}
	void removeAllRecord(){
	  records.clear();
	}

	bool load(const string file_name);
	bool save(const string file_name)const;
	
  private:
	vector<OneDragRecord> records;
	
  };
  
  typedef boost::shared_ptr<DragTrajectoryRecord> pDragTrajectoryRecord;
  
}//end of namespace

#endif /*_DRAGTRAJECTORYRECORD_H_*/
