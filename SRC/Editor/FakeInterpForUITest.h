#ifndef _FAKEINTERPFORUITEST_H_
#define _FAKEINTERPFORUITEST_H_

#include <boost/shared_ptr.hpp>
#include <BaseInterpolator.h>
#include <assertext.h>
#include <JsonFilePaser.h>
#include <TetMesh.h>
#include <MatrixIO.h>
#include <MatrixTools.h>
#include <Log.h>
using namespace UTILITY;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class FakeInterpForUITest this is a fake interpolator, which is just used
   * for testing the UI.
   * 
   */
  class FakeInterpForUITest: public BaseInterpolator{
	
  public:
	FakeInterpForUITest(){
	  WARN_LOG("the fake interpolator for UI testing will be used.");
	}
	bool init (const string init_filename){

	  JsonFilePaser json_f;
	  bool succ = json_f.open(init_filename);

	  TetMesh tet_mesh;
	  { // load vol mesh
		string vol_filename;
		succ = json_f.readFilePath("vol_filename",vol_filename);
		if (succ) succ = tet_mesh.load(vol_filename);
	  }

	  if (succ){
		// read input sequence
		int zero_input_animation = 0;
		if (json_f.read("zero_input_animation",zero_input_animation) && zero_input_animation > 0){
		  _inputU.resize(zero_input_animation);
		  for (size_t i = 0; i < _inputU.size(); ++i){
			_inputU[i].resize(tet_mesh.nodes().size()*3);
			_inputU[i].setZero();
		  }
		}else{
		  	succ = false;
			string fullinput;
			if(json_f.readFilePath("full_input_animation",fullinput)){
			  MatrixXd Uin;
			  succ = EIGEN3EXT::load(fullinput,Uin);
			  EIGEN3EXT::convert(Uin,_inputU);
			}
		}
	  }
	  return succ;
	}
	bool editable(const int frame_id)const{
	  return true;
	}
	void setConGroups(const int frame_id,const vector<set<int> >&group, 
							  const Eigen::Matrix<double,3,-1> &pc){}
	void setUc(const int frame_id, const Eigen::Matrix<double,3,-1> &pc){}
	void removeAllPosCon(){}
	void removePartialCon(const int frame_id){}
	bool interpolate (){
	  return true;
	}
	const VectorXd& getInterpU(const int frame_id){
	  assert_in(frame_id,0,_inputU.size()-1);
	  return _inputU[frame_id];
	}
	const VectorXd& getInputU(const int frame_id){
	  assert_in(frame_id,0,_inputU.size()-1);
	  return _inputU[frame_id];
	}
	const VectorXd& getReducedEdits(const int frame_id)const{
	  assert_in(frame_id,0,_inputU.size()-1);
	  return _inputU[frame_id];
	}
	int getT()const{
	  return _inputU.size();
	}
	int reducedDim()const{
	  return _inputU.size()>0? _inputU[0].size():0;
	}
	string getName()const{
	  return "fake interpolator for UI testing";
	}

  private:
	vector<VectorXd> _inputU; 
  };
  
  typedef boost::shared_ptr<FakeInterpForUITest> pFakeInterpForUITest;
  
}//end of namespace

#endif /*_FAKEINTERPFORUITEST_H_*/
