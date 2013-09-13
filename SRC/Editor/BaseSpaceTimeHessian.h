#ifndef _BASESPACETIMEHESSIAN_H_
#define _BASESPACETIMEHESSIAN_H_

#include <vector>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <assertext.h>
#include <boost/shared_ptr.hpp>
using namespace Eigen;
using namespace std;

namespace IEDS{

  typedef Eigen::Triplet<double> E_Triplet;
  
  /**
   * @class BaseSpaceTimeHessian base class for calculating the hessian of the
   * spacetime function.
   * 
   */
  class BaseSpaceTimeHessian{
	
  public:
	BaseSpaceTimeHessian(){
	  T = 3;
	  h = 0.0f;
	  alpha_k = 0.0f;
	  alpha_m = 0.0f;
	}

	// initialize
	void setTotalFrame(const int T){
	  assert_ge(T,3); // 3 frames at least.
	  this->T = T;
	}
	void setTimeStep(const double h){
	  assert_gt(h,0.0f);
	  this->h = h;
	}
	void setStiffnessDamping(const double s){
	  assert_ge(s,0.0f);
	  alpha_k = s;
	}
 	void setMassDamping(const double s){
	  assert_ge(s,0.0f);
	  alpha_m = s;
	}
	void setEigenValues(const VectorXd &eigen_values){
	  this->eigen_values = eigen_values;
	}

	int reducedDim()const{
	  return eigen_values.size();
	}
	int totalFrames()const{
	  return T;
	}
	int totalDim()const{
	  return totalFrames()*reducedDim();
	}

	bool compute_H(SparseMatrix<double> &H){

	  vector<E_Triplet> H_triplet;
	  const bool succ = this->compute_HTriplet(H_triplet);
	  if (succ){
		H.resize(totalDim(),totalDim());
		H.reserve(H_triplet.size());
		H.setFromTriplets(H_triplet.begin(), H_triplet.end());
	  }
	  return succ;
	}
	// compute hessian matrix, should be implemented in the subclass.
	virtual bool compute_HTriplet(vector<E_Triplet> &H_triplet) = 0;

  protected:
	int T; // total frames
	double h; //time step
	double alpha_k; // stiffness damping
	double alpha_m; // mass damping	
	VectorXd eigen_values;

  };
  
  typedef boost::shared_ptr<BaseSpaceTimeHessian> pBaseSpaceTimeHessian;
  
}//end of namespace

#endif /*_BASESPACETIMEHESSIAN_H_*/
