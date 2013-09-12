#ifndef _NODEWARPER_H_
#define _NODEWARPER_H_

#include <assertext.h>
#include <vector>
#include <eigen3/Eigen/Dense>
typedef std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > VectorV3;

namespace LSW_WARPING{
  
  /**
   * @class NodeWarper construct the well deformed shape by using the
   * per-node-warping approach introduced by Modal-Warping[TVCG2005].
   * 
   */
  class NodeWarper{
	
  public:
	static void warp(const VectorV3 &w, const VectorXd &p, VectorXd &u){
	  
	  assert_eq((int)w.size()*3,p.size());
	  u.resize(p.size());
	  Vector3d u3;
	  for (size_t i = 0; i < w.size(); ++i){
		const Vector3d p3 = p.segment(i*3,3);
		warp(w[i],p3,u3);
		u.segment(i*3,3) = u3;
	  }
	}

	static void warp(const Vector3d &w, const Vector3d &p, Vector3d &u){
	  
	  if (w.norm() <= 1e-8){/// @todo
		u = p;
	  }else{
		const Matrix3d _R = warpMat(w);
		u = _R*p;
	  }
	}

	static void warpMat(const VectorV3 &w, SparseMatrix<double> &warpedR){

	  vector<Matrix3d> Rs(w.size());
	  for (size_t i = 0; i < w.size(); ++i)
		Rs[i] =  warpMat(w[i]);
	  assembleR3x3(Rs, warpedR);
	}

	static void rotMat(const VectorV3&w,SparseMatrix<double> &R){
	  
	  vector<Matrix3d> Rs(w.size());
	  for (size_t i = 0; i < w.size(); ++i)
		Rs[i] =  rotMat(w[i]);
	  assembleR3x3(Rs, R);
	}

	static Matrix3d warpMat(const Vector3d &w){

	  const Matrix3d I = Matrix3d::Identity(3,3);
	  const double nw = w.norm();
	  if (nw < 1e-8){
		return I;
	  }else{
		const Vector3d _w = w.normalized();
		const Matrix3d wx = SkewMat(_w);
		return I + wx*( (1.0f-cos(nw))/nw ) + wx*wx*( 1.0f-sin(nw)/nw );
	  }
	}

	static Matrix3d rotMat(const Vector3d &w){
	  
	  const Matrix3d I = Matrix3d::Identity(3,3);
	  const double nw = w.norm();
	  if (nw < 1e-8){
		return I;
	  }else{
		const Vector3d _w = w.normalized();
		const Matrix3d wx = SkewMat(_w);
		return I + wx*sin(nw) + wx*(wx*( 1.0f-cos(nw) ));
	  }
	}

	static Matrix3d SkewMat(const Vector3d &_w){

	  Matrix3d wx;
	  wx.setZero();
	  
	  wx(0,1) = -_w[2];
	  wx(0,2) = _w[1];
	  wx(1,2) = -_w[0];

	  wx(1,0) = -wx(0,1);
	  wx(2,0) = -wx(0,2);
	  wx(2,1) = -wx(1,2);

	  return wx;
	}

	static void assembleR3x3(const vector<Matrix3d> &Rs, SparseMatrix<double> &M){
	  
	  std::vector<Eigen::Triplet<double> > triplet;
	  triplet.reserve(Rs.size()*9);

	  for (size_t i = 0; i < Rs.size(); ++i){

		const Matrix3d &R = Rs[i];
		for (int r = 0; r < 3; ++r){
		  for (int c = 0; c < 3; ++c){
			triplet.push_back( Eigen::Triplet<double>(i*3+r,i*3+c,R(r,c)) );
		  }
		}
	  }
	  
	  M.resize (Rs.size()*3,Rs.size()*3);
	  M.reserve (triplet.size());
	  M.setFromTriplets( triplet.begin(),triplet.end());
	}
	
  };
  
}//end of namespace

#endif /*_NODEWARPER_H_*/
