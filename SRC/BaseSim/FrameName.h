#ifndef _FRAMENAME_H_
#define _FRAMENAME_H_

#include <string>
using namespace std;

namespace LSW_RENDERMESH{
  
  /**
   * @class FrameName store the name of a frame of an animation sequence. The
   * name of each frame always has a number which is indicate the frame number
   * of this file, such as "beam99.obj" (frame 99 of the beam animation). It is
   * also provide a method to compare between two file, which is used to sort
   * the file sequence.
   * 
   */
  class FrameName{
	
  public:
	FrameName(){
	  frame_num = -1;
	}
	FrameName(const string &_filename){
	  setFileName(_filename);
	}
	void setFileName(const string &_filename);
	bool operator < (const FrameName& f)const;
	int getFrameNum()const{
	  return frame_num;
	}
	string toString()const{
	  return filename;
	}
	
  protected:
	int grepNumber(const string &filename)const;
	
  private:
	string filename;
	int frame_num;
  };
  
}//end of namespace

#endif /*_FRAMENAME_H_*/
