#ifndef _FASTANIEDITOR_H_
#define _FASTANIEDITOR_H_

#include <string>
#include <vector>
using namespace std;
#include <BaseSpaceTimeHessian.h>

namespace IEDS{
  
  typedef vector<Eigen::Triplet<double> > VecT;

  /**
   * @class FastAniEditor the main entry of the fast animation editing system,
   * providing a set of functions to initialize the editing system and
   * manipulate the sequence, producing the edited results.
   * 
   * @see ReadMe
   */
  class FastAniEditor{
	
  public:
	FastAniEditor();

	/// read these data from the init file: h, alpha_k, alpha_m, r,T,eigen_values
	bool initialize(const string &ini_file);
	void setBoundayFrames(const vector<int> &boundary_frames);
	void setTotalFrame(const int T);
	void setTimeStep(const double h);
	void setStiffnessDamping(const double s);
	void setMassDamping(const double s);
	void setEigenValues(const VectorXd &eigen_values);

	// calculate Hessian Matrix
	bool precompute();

	// manipulate
	bool editable(const int frame_id)const;// check if the frame can be edited.
	bool isBoundaryFrame(const int frame_id)const;
	void setPosCon(const vector<int> &con_frame_id,
				   const vector<MatrixXd> &C,
				   const vector<VectorXd> &delta_uc);
	bool isRemoved(const int frame_id)const;

	// apply edits
	bool applyEdits();

	// get result
	const VectorXd &getEditResult()const{
	  return Z;
	}

	const SparseMatrix<double> &get_A()const{
	  return A;
	}
	int reducedDim()const;
	int numBoundaryFrames()const;
	int totalFrames()const;
	int totalDim()const;
	const VectorXd &getEigenvalues()const{
	  assert(hessian_comp);
	  return hessian_comp->getEigenvalues();
	}

  protected:
	int removedFrameNumBefore(const int frame_id)const;
	void correctBoundaryFrames();
	
  private:
	pBaseSpaceTimeHessian hessian_comp;
	// tripletlist of the hessian matrix H of the energy function.
	VecT H_triplet;

	/**
	 * A X = P, where X = [Z , lambda]^t, and P = [O , P^m]^t , and lambda is
	 * the Lagrange multiplier, and P^m is the postion constraints, and A is
	 * | H C^t |
	 * | C  O  |. 
	 * Both A and P are initialized by calling setPosCon(...).
	 * @see Barbic's siggraph 2012 paper.
	 * @note A is a symetric matrix and is full stored.
	 * @note as frames in boundary_frames is set as boudary conditions, and won't
	 * be changed, the corresponding rows and columns of A will be removed.
	 */
	SparseMatrix<double> A; 
	VectorXd P;
	// boundary conditions, will be removed from A, see BSG2012. By default,
	// frame 0,1 will be set as boundary conditions.
	vector<int> boundary_frames; 

	/// the increamental displacements in the MA coordinates, the length of
	/// which is (T-1)*r.
	VectorXd Z; 
  };
  
  typedef boost::shared_ptr<FastAniEditor> pFastAniEditor;
  
}//end of namespace

#endif /*_FASTANIEDITOR_H_*/
