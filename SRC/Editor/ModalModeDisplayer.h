#ifndef _MODALMODEDISPLAYER_H_
#define _MODALMODEDISPLAYER_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include <Log.h>
#include <JsonFilePaser.h>
using namespace Eigen;
using namespace std;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class ModalModeDisplayer display modal modes.
   * 
   */
  class ModalModeDisplayer{
	
  public:
	ModalModeDisplayer(){
	  show_modal_modes = false;
	  modal_scale = 5000.0f;
	  number_modes_to_show = 20;
	  displayed_mode = -1;
	}
	bool isShowModalModes()const{
	  return show_modal_modes;
	}
	bool initialize(const string init_filename){
	  
	  UTILITY::JsonFilePaser json_f;
	  bool succ = json_f.open(init_filename);
	  if (succ){
		if (json_f.read("show_modal_modes",show_modal_modes) && show_modal_modes){
		  WARN_LOG("modal modes will be shown instead of interpoaltion results.");
		  if(!json_f.read("modal_scale",modal_scale)){
			modal_scale = 5000.0f;
		  }
		  if(!json_f.read("number_modes_to_show",number_modes_to_show)){
			number_modes_to_show = 20;
		  }
		}
	  }else{
		show_modal_modes = false;
	  }
	  return succ;
	}
	bool showModalModes(const VectorXd &eigenvalues, vector<VectorXd> &z){

	  if (!show_modal_modes){
		return false;
	  }
	  
	  const int T = z.size();
	  if (T > 0){
		const int r = z[0].size();
		assert_gt(r,0);
		const double s = modal_scale;
		if ( displayed_mode>=number_modes_to_show || displayed_mode >= r ){
		  displayed_mode = -1;
		}
		displayed_mode ++;
		const double zT = s/sqrt(eigenvalues[displayed_mode]);
		cout << "mode = " << displayed_mode << ", z[T] = "<< zT << endl;
		for (size_t f = 0; f < T; ++f){
		  z[f].setZero();
		  z[f][displayed_mode]=zT*f/T;
		}
	  }

	  return true;
	}
	
  private:
	bool show_modal_modes;
	double modal_scale;
	int number_modes_to_show;
	int displayed_mode;
  };

}//end of namespace

#endif /*_MODALMODEDISPLAYER_H_*/
