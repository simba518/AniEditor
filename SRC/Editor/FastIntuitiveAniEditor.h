#ifndef _FASTINTUITIVEANIEDITOR_H_
#define _FASTINTUITIVEANIEDITOR_H_

#include <string>
#include <vector>
#include <BaseSpaceTimeHessian.h>
#include <SpaceTimeHessianInverse.h>
using namespace std;

namespace IEDS{
  
  /**
   * @class FastIntuitiveAniEditor implement a fast intuitive animation editing
   * method. The algorithm is similar to FastAniEditor, but with much efficient
   * than FastAniEditor.
   */
  class FastIntuitiveAniEditor{
	
  public:
	// init
	FastIntuitiveAniEditor(){
	  boundary_frames.push_back(0);
	  boundary_frames.push_back(1);
	  T = 0;
	}
	bool initialize(const string &ini_file);
	void setBoundayFrames(const vector<int> &b_frames){
	  this->boundary_frames = b_frames;
	}
	void setTotalFrame(const int T){
	  assert_ge(T,2);
	  this->T = T;
	  correctBoundaryFrames();
	}

	// assemble H and calculate inv(H).
	void precompute();
	
	// manipulate
	bool isBoundaryFrame(const int frame_id)const{
	  for (size_t i = 0; i < boundary_frames.size(); ++i){
		if (frame_id == boundary_frames[i]){
		  return true;
		}
	  }
	  return false;
	}
	void setConstrainedFrames(const vector<int>&con_frame_id);
	void clearConstraints();
	void computeConstrainedFrames(const vector<MatrixXd> &Cs, 
								  const VectorXd &uc,
								  VectorXd &zc);
	void updateZ();

	// get result
	bool hasConstraints()const{
	  return invHcSet.size() > 0;
	}
	const VectorXd &getEditResult()const{
	  return z;
	}
	int reducedDim()const{
	  return eigen_values.size();
	}
	int totalFrame()const{
	  return T;
	}
	int numBoundaryFrames()const{
	  return (int)boundary_frames.size();
	}
	const VectorXd &getEigenvalues()const{
	  return eigen_values;
	}
	
  protected:
	int removedFrameNumBefore(const int frame_id)const;
	void correctBoundaryFrames();

  private:
	vector<int> boundary_frames; // boundary frames, will be removed from H.
	SparseMatrix<double> invH; // inverse of H, with boundary frames removed.
	SparseMatrix<double> invHc;
	vector<SparseMatrix<double> > invHcSet;
	SparseMatrix<double> invHcZ;
	VectorXd C_lambda; // C*(-lambda)
	VectorXd z; // total edits.
	int T; // total frames
	double h; //time step
	double alpha_k; // stiffness damping
	double alpha_m; // mass damping	
	VectorXd eigen_values;
  };

}//end of namespace

#endif /*_FASTINTUITIVEANIEDITOR_H_*/
