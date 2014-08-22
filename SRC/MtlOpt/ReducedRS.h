#ifndef _REDUCEDRS_H_
#define _REDUCEDRS_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <ComputeBj.h>
#include <MatrixIO.h>
#include <RSCoordComp.h>
using namespace std;
using namespace Eigen;
using namespace LSW_WARPING;

namespace MTLOPT{
  
  /**
   * @class ReducedRS Efficient reduced RS method.
   * 
   */
  class ReducedRS{
	
  public:
	ReducedRS(){}
	inline void checkDimensions(){
	  assert_eq(B.rows()%3,0);
	  assert_eq(B.cols(),P.rows());
	  assert_eq(P.cols(),refY.rows());
	  assert_eq(refY.rows()%9,0);
	  assert_eq(hatW.rows(),P.cols());
	  assert_eq(SqrtV.size()*9,P.cols());
	}
	inline void warp(const VectorXd &z,const int f, VectorXd &u){

	  checkDimensions();
	  assert_in(f,0,refY.cols()-1);
	  static VectorXd b;
	  const VectorXd y = hatW*z + refY.col(f);
	  b.resize(y.size());
	  computeB(y,f,b);
	  assert_eq(P.cols(),b.size());
	  const VectorXd q = P*b;
	  u = B*q;
	}

	void setBP(const MatrixXd &b,const MatrixXd &p){
	  B = b;
	  P = p;
	}
	void setHatW(const MatrixXd &hw){
	  hatW = hw;
	}
	template<class VECTOR>
	void setSqrtV(const VECTOR &v, const bool devide_9){
	  if (devide_9){
		assert_eq(v.size()%9,0);
		SqrtV.resize(v.size()/9);
		for (int i = 0; i < v.size(); i+=9)
		  SqrtV[i/9] = v[i];
	  }else{
		SqrtV.resize(v.size());
		for (int i = 0; i < v.size(); ++i)
		  SqrtV[i] = v[i];
	  }
	}
	void setRefY(const MatrixXd &y){
	  refY = y;
	}
	void setRefY(const vector<VectorXd> &y){
	  assert_gt(y.size(),0);
	  refY.resize(y[0].size(),y.size());
	  for (size_t i = 0; i < y.size(); ++i){
		assert_eq_ext(y[i].size(),refY.rows(),"i="<<i);
		refY.col(i) = y[i];
	  }
	}

	const MatrixXd &getB()const{
	  return B;
	}
	const MatrixXd &getP()const{
	  return P;
	}
	const MatrixXd &getHatW()const{
	  return hatW;
	}
	const VectorXd &getSqrtV()const{
	  return SqrtV;
	}
	const MatrixXd &getRefY()const{
	  return refY;
	}

	inline void computeB(const VectorXd &yi,int f, VectorXd &b)const{
	  // /// @todo optimize using multi-threads.
	  VectorXd y(yi.size());
	  y.setConstant(1e-14);
	  y += yi;
	  for (int i = 0; i < y.size(); i+=9){
		ComputeBj::compute(&y[i],SqrtV[i/9],&b[i]);
	  }
	}
	bool write(const string filename)const{
	  bool succ = EIGEN3EXT::write(filename+"_B.b", B);
	  succ &= EIGEN3EXT::write(filename+"_hatW.b", hatW);
	  succ &= EIGEN3EXT::write(filename+"_P.b", P);
	  succ &= EIGEN3EXT::write(filename+"_refY.b", refY);
	  succ &= EIGEN3EXT::write(filename+"_SqrtV.b", SqrtV);
	  return succ;
	}
	bool load(const string filename){
  	  bool succ = EIGEN3EXT::load(filename+"_B.b", B);
	  succ &= EIGEN3EXT::load(filename+"_hatW.b", hatW);
	  succ &= EIGEN3EXT::load(filename+"_P.b", P);
	  succ &= EIGEN3EXT::load(filename+"_refY.b", refY);
	  succ &= EIGEN3EXT::load(filename+"_SqrtV.b", SqrtV);
	  return succ;
	}

  private:
	// cubature applied for all data bellow.
	MatrixXd B; 
	MatrixXd hatW;
	MatrixXd P; // (B^t*A*B)^{-1}B^T*((VG)^t)
	MatrixXd refY;
	VectorXd SqrtV;
  };
  
  typedef boost::shared_ptr<ReducedRS> pReducedRS;
  typedef boost::shared_ptr<const ReducedRS> pReducedRS_const;

  // convert a full u to z according the reduced rs method.
  // @bug no refY used yet, e.g can only be used for static input sequence.
  class ReducedRSUnwarp{
	
  public:
	void setHatW(const MatrixXd &hw){
	  hatW = hw;
	}
	void setSubG(const SparseMatrix<double> &sG){
	  subG = sG;
	}
	void convert(const MatrixXd &U, MatrixXd &Z)const{
	  checkDimensions();
	  MatrixXd Y;
	  RSCoordComp::constructWithoutWarp(subG,U,Y);
	  Z = (hatW.transpose()*hatW).ldlt().solve(hatW.transpose()*Y);
	}
	void checkDimensions()const{
	  assert_gt(hatW.cols(),0);
	  assert_gt(hatW.rows(),0);
	  assert_eq(subG.rows(),hatW.rows());
	}
	
  private:
	MatrixXd hatW;
	SparseMatrix<double> subG;
  };
  
}//end of namespace

#endif /*_REDUCEDRS_H_*/
