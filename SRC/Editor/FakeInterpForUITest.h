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

	  TRACE_FUN();
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

	  _group.resize(getT());
	  _pc.resize(getT());
	  _ouputU = _inputU;
	  return succ;
	}
	bool editable(const int frame_id)const{
	  assert_in(frame_id,0,getT()-1);
	  return true;
	}
	void setConGroups(const int frame_id,const vector<set<int> >&group, 
							  const Eigen::Matrix<double,3,-1> &pc){
	  assert_in(frame_id,0,getT());
	  int nodes = 0;
	  for (int i = 0; i < group.size(); ++i){
		assert_gt(group[i].size(),0);
		_group[frame_id] = group;
		_pc[frame_id] = pc;
		nodes += group[i].size();
	  }
	  assert_eq(pc.cols(),nodes);
	}
	void setUc(const int frame_id, const Eigen::Matrix<double,3,-1> &pc){
	  assert_in(frame_id,0,getT()-1);
	  assert_gt(pc.cols(),0);
	  assert_eq(_pc[frame_id].cols(),pc.cols());
	  _pc[frame_id] = pc;
	}
	void removeAllPosCon(){
	  BOOST_FOREACH(vector<set<int> > &g, _group){
		g.clear();
	  }
	  for (size_t i = 0; i < getT(); ++i){
		_pc[i].resize(3,0);
	  }
	}
	void removePartialCon(const int frame_id){
	  assert_in(frame_id,0,getT()-1);
	  _group[frame_id].clear();
	  _pc[frame_id].resize(3,0);
	}
	bool interpolate (){

	  for (int f = 0; f < getT(); ++f){
		const vector<set<int> > &group = _group[f];
		const Eigen::Matrix<double,3,-1> &pc = _pc[f];
		int nodes = 0;
		for (int j = 0; j < group.size(); ++j){
		  BOOST_FOREACH(const int i, group[j]){
		    assert_in(i*3,0,_inputU[f].size()-2);
			assert_lt(nodes,pc.cols());
		    _ouputU[f].segment<3>(i*3) = pc.col(nodes);
			nodes++;
		  }
		}
	  }
	  return true;
	}
	const VectorXd& getInterpU(const int frame_id){
	  assert_in(frame_id,0,_ouputU.size()-1);
	  return _ouputU[frame_id];
	}
	const VectorXd& getInputU(const int frame_id){
	  assert_in(frame_id,0,_inputU.size()-1);
	  return _inputU[frame_id];
	}
	const VectorXd& getReducedEdits(const int frame_id)const{
	  assert_in(frame_id,0,_ouputU.size()-1);
	  return _ouputU[frame_id];
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
	vector<VectorXd> _ouputU;
	vector<vector<set<int> > > _group;
	vector<Eigen::Matrix<double,3,-1> > _pc;
  };
  
  typedef boost::shared_ptr<FakeInterpForUITest> pFakeInterpForUITest;
  
}//end of namespace

#endif /*_FAKEINTERPFORUITEST_H_*/
