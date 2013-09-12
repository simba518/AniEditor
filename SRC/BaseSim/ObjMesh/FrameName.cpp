#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>
#include <assertext.h>
#include <Log.h>
#include "FrameName.h"
using namespace boost;
using namespace LSW_RENDERMESH;

#define IS_DIGIT(one_char)						\
  (one_char <= '9' && one_char >= '0')

#define IS_NOT_DIGIT(one_char)					\
  (!(IS_DIGIT(one_char)))

bool FrameName::operator < (const FrameName& f)const{
  
  return (this->getFrameNum() < f.getFrameNum());
}

void FrameName::setFileName(const string &_filename){

  this->filename = _filename;
  this->frame_num = grepNumber(filename);
}

int FrameName::grepNumber(const string &filename)const{
  
  // file name without directory
  const string fname = filesystem::path(filename).filename().string(); 
  int number_start = -1;
  for (int i = 0; i < (int)fname.size(); ++i){
    if (IS_DIGIT(fname[i])){
	  number_start = i;
	  break;
	}
  }

  if (number_start < 0){
	ERROR_LOG("there is no number in the filename: " << filename);
	return -1;
  }

  int number_end = number_start+1;
  for (int i = number_start+1; i < (int)fname.size(); ++i){
	if (IS_NOT_DIGIT(fname[i])){
	  number_end = i;
	  break;
	}
  }

  const int len = number_end-number_start;
  const int number = lexical_cast<int>(fname.substr(number_start,len));

  return number;
}
