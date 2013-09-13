#include <fstream>
#include "DragTrajectoryRecord.h"
using namespace std;
using namespace LSW_ANI_EDITOR;

bool DragTrajectoryRecord::load(const string file_name){

  this->removeAllRecord();
  bool succ = true;
  ifstream in(file_name.c_str());
  if (in.is_open()){
	int total_record = 0;
	in >> total_record;
	assert_ge (total_record,0);
	for (int i = 0; i < total_record; ++i){

	  // laod frame number
	  int frame_num = 0;
	  in >> frame_num;

	  // load uc
	  int uc_len = 0;
	  in >> uc_len;
	  assert_ge (uc_len, 0);
	  VectorXd uc(uc_len);
	  for (int j = 0; j < uc_len; ++j){
		in >> uc[j];
	  }

	  // record drag trajectory.
	  this->record(frame_num, uc);
	}
  }else{
	succ = false;
  }
  return succ;  
}

bool DragTrajectoryRecord::save(const string file_name)const{

  bool succ = true;
  ofstream out(file_name.c_str());
  if (out.is_open()){
	
	out << totalRecord() << endl;
	for (int i = 0; i < totalRecord(); ++i){

	  const OneDragRecord &one_record = getRecord(i);
	  const int frame_num = one_record.getFrameNum();
	  const VectorXd &uc = one_record.getUc();
	  out << frame_num << endl;
	  out << uc.size() << endl;
	  for (int j = 0; j < uc.size(); ++j){
		out << uc[j] << "\t";
	  }
	}
  }else{
	succ = false;
  }
  return succ;
}
