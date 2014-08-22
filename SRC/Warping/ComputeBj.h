#ifndef _COMPUTEBJ_H_
#define _COMPUTEBJ_H_

#include <memory.h>
#include <iostream>
#include <eigen3/Eigen/Dense>
#include <boost/shared_ptr.hpp>
#include <assertext.h>
using namespace Eigen;
using namespace std;

namespace LSW_WARPING{

  typedef Eigen::Matrix<double, 3, 3, Eigen::RowMajor> mat3r;
  
  /**
   * @class ComputeBj compute bj = sqrt(Vj)*(exp(yj_w)*(I+S(yj_s)) - I).
   * @see Barbic's siggraph 2012.
   * @todo optimize
   */
  class ComputeBj{
	
  public: 
	// compute bj
	/// @note yj and bj are 9x1 vectors expressed as a 3x3 row-major matrix.
	static void compute(const double *yj,const double Sqrt_Vj, double *bj){
  
	  assert_gt(Sqrt_Vj,0);

	  // the symtric part
	  mat3r SI;
	  compute_S(yj, SI);
	  SI(0,0) += 1.0f;
	  SI(1,1) += 1.0f;
	  SI(2,2) += 1.0f;
  
	  // the antisymmetric part
	  mat3r Exp_W;
	  compute_exp(yj, Exp_W);

	  // bj
	  Eigen::Matrix<double, 3, 3, RowMajor> m;
	  m = Sqrt_Vj*(Exp_W*SI);
	  m(0,0) -= Sqrt_Vj;
	  m(1,1) -= Sqrt_Vj;
	  m(2,2) -= Sqrt_Vj;
	  memcpy(bj,m.data(),sizeof(double)*9);
	}
	static void compute(const VectorXd &yj,const double Sqrt_Vj,VectorXd &bj){
	  assert_eq(yj.size(),9);
	  bj.resize(9);
	  compute(&yj[0],Sqrt_Vj,&bj[0]);
	}

	// compute only the symetric part: S(yj_s)
	static mat3r &compute_S(const double *yj, mat3r &SI){

	  SI(0,0) = yj[0+3];
	  SI(1,0) = yj[1+3];
	  SI(2,0) = yj[2+3];
	  SI(1,1) = yj[3+3];
	  SI(2,1) = yj[4+3];
	  SI(2,2) = yj[5+3];

	  SI(0,1) = SI(1,0);
	  SI(0,2) = SI(2,0);
	  SI(1,2) = SI(2,1);
	  return SI;
	}
	static mat3r &compute_S(const VectorXd &yj, mat3r &m){
	  assert_ge (yj.size(), 3);
	  return compute_S(&yj[0],m);
	}
	static mat3r compute_S(const double *yj){
	  mat3r m;
	  return compute_S(yj, m);
	}
	static mat3r compute_S(const VectorXd &yj){
	  mat3r m;
	  return compute_S(yj, m);
	}

	// compute only the rotational part: exp(yj_w)
	// @todo optimize using maxima
	static mat3r &compute_exp(const double *yj, mat3r &Exp_W){

	  const double theta = sqrt(yj[0]*yj[0] + yj[1]*yj[1] + yj[2]*yj[2]);
	  assert_ge(theta,0.0f);

	  /// @bug remove the if condition.
	  if(theta < 1e-8){ // the rotation angle is too small
		Exp_W.setZero();
		Exp_W(0,0) = 1;  Exp_W(1,1) = 1;  Exp_W(2,2) = 1;
	  }else{
		const double x = yj[0]/theta;
		const double y = yj[1]/theta;
		const double z = yj[2]/theta;
		const double c = 1.0f - cos(theta);
		Exp_W(0,0) = cos(theta)+x*x*c ;      Exp_W(0,1) = x*y*c-z*sin(theta);     Exp_W(0,2) = x*z*c + y* sin(theta);  // row 0
		Exp_W(1,0) = x*y *c + z*sin(theta) ; Exp_W(1,1) = cos(theta) + y*y *c ;   Exp_W(1,2) =  y* z*c - x*sin(theta); // row 1
		Exp_W(2,0) = z*x*c - y* sin(theta) ; Exp_W(2,1) = z*y *c + x*sin(theta) ; Exp_W(2,2) =  cos(theta) + z*z *c;   // row 2
	  }
	  return Exp_W;
	}
	static mat3r &compute_exp(const VectorXd &yj, mat3r &m){
	  assert_ge (yj.size(), 3);
	  return compute_exp(&yj[0],m);
	}
	static mat3r compute_exp(const double *yj){
	  mat3r m;
	  return compute_exp(yj, m);
	}
	static mat3r compute_exp(const VectorXd &yj){
	  mat3r m;
	  return compute_exp(yj, m);
	}
	
  };
  
  typedef boost::shared_ptr<ComputeBj> pComputeBj;
  
}//end of namespace

#endif /*_COMPUTEBJ_H_*/
