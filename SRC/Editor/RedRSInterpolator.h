#ifndef _REDRSINTERPOLATOR_H_
#define _REDRSINTERPOLATOR_H_

#include <FastIntuitiveAniEditor.h>
#include <WarpInterpolator.h>
using namespace IEDS;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class RedRSInterpolator intuitive animation interpolator based on reduced
   * RS method.
   * 
   */
  class RedRSInterpolator: public WarpInterpolator{
	
  public:
	bool init (const string init_filename);
	bool editable(const int frame_id)const{
	  return !(anieditor.isBoundaryFrame(frame_id));
	}
	void setConGroups(const int f_id,const vector<set<int> >&group,const Matrix<double,3,-1> &uc){
	  vector<set<int > > group_rlst;
	  VectorXd uc_rlst;
	  splitAllConstraints(group,uc,group_rlst,uc_rlst);
	  addConGroups(f_id,group_rlst,uc_rlst);
	  anieditor.setConstrainedFrames(con_frame_id);
	}
	bool interpolate ();
	string getName()const{
	  return string("Reduced RS");
	}
	void print(const string data_name)const{
	  if(string("eigenvalues") == data_name){
		cout<< "eigenvalues: \n" << anieditor.getEigenvalues() << endl;
	  }else{
		BaseInterpolator::print(data_name);
	  }
	}

  private:
	FastIntuitiveAniEditor anieditor;
  };
  
  typedef boost::shared_ptr<RedRSInterpolator> pRedRSInterpolator;
  
}//end of namespace

#endif /* _REDRSINTERPOLATOR_H_*/
