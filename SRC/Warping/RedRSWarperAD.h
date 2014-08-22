#ifndef _REDRSWARPERAD_H_
#define _REDRSWARPERAD_H_

#include <set>
#include <boost/shared_ptr.hpp>
#include <symbolic/fx/sx_function.hpp>
#include <algorithm>
#include <SparseMatrixTools.h>
#include <RSCoordComp.h>
#include <ComputeBjAD.h>
#include <CASADITools.h>
#include <MapMA2RS.h>
#include <RS2Euler.h>
using namespace std;
using CasADi::SX;
using CasADi::SXMatrix;

namespace LSW_WARPING{

  /**
   * @class RedRSWarperGrad compute q and u of the reduced RS method.
   * 
   */
  class RedRSWarperAD{
	
  public:
	// initialize
	RedRSWarperAD(){
	  
	}
	RedRSWarperAD(const RS2Euler &rs2euler,const MatrixXd &NLBasis,const MatrixXd &LB){
	  init(rs2euler,NLBasis,LB);
	  convert2Symbols();
	}
	RedRSWarperAD(const RS2Euler &rs2euler,const MatrixXd &NLBasis,const MatrixXd &LBasis,
				  const vector<int> &selPoints, const vector<double> &cubWeights){
	  init(rs2euler,NLBasis,LBasis);
	  setCubature(selPoints, cubWeights);
	  convert2Symbols();
	}

	// set directly
	void setBP(const MatrixXd &b,const MatrixXd &p){
	  B = b;
	  P = p;
	  CASADI::convert(B,Bs);
	  CASADI::convert(P,Ps);
	}
	void setHatW(const MatrixXd &hw){
	  hatW = hw;
	  CASADI::convert(hatW,hatWs);
	}
	template<class VECTOR>
	void setSqrtV(const VECTOR &v){

	  Sqrt_V.clear();
	  Sqrt_Vs.clear();
	  Sqrt_V.reserve(v.size());
	  Sqrt_Vs.reserve(v.size());
	  for (int i = 0; i < v.size(); ++i){
		Sqrt_V.push_back(v[i]);
		Sqrt_Vs.push_back(v[i]);
	  }
	}
	void checkDimensions(){
	  assert_eq(B.rows()%3,0);
	  assert_eq(B.cols(),P.rows());
	  assert_eq(hatW.rows(),P.cols());
	  assert_eq(Sqrt_V.size()*9,P.cols());
	  assert_eq(Sqrt_Vs.size()*9,P.cols());
	}

	// warp
	void warp(const VectorXd &z,VectorXd &u){
	  const VectorXd input_y = VectorXd::Zero(hatW.rows());
	  warp(input_y,z,u);
	}
	void warp(const VectorXd &z,VectorXd &u,const vector<int> &nodes){
	  const SXMatrix zs = CASADI::convert(z);
	  SXMatrix us;
	  warp(zs,us,nodes);
	  u = CASADI::convert2Vec<double>(us);
	  assert_eq(u.size(),nodes.size()*3);
	}
	void warp(const SXMatrix &z,SXMatrix &u,const vector<int> &nodes){
	  const VectorXd y = VectorXd::Zero(hatWs.size1());
	  warp(y,z,u,nodes);
	}
	void warp(const VectorXd &input_y,const SXMatrix &z,SXMatrix &u,const vector<int> &nodes){

	  SXMatrix q;
	  computeQ(input_y,z,q);
	  set<int> s;
	  for (int i = 0; i < nodes.size(); ++i) s.insert(nodes[i]);
	  Eigen::SparseMatrix<double> P;
	  EIGEN3EXT::genReshapeMatrix(B.rows(),3,s,P,false);
	  const SXMatrix subB = CASADI::convert((MatrixXd)(P*B));
	  assert_eq(subB.size1(),nodes.size()*3);
	  assert_eq(subB.size2(),q.size1());
	  u = subB.mul(q);
	}
	void warp(const VectorXd &input_y,const SXMatrix &z,SXMatrix &u){
	  SXMatrix q;
	  computeQ(input_y,z,q);
	  u = Bs.mul(q);
	}
	void warp(const VectorXd &input_y,const VectorXd &z,VectorXd &u){
	  VectorXd q;
	  computeQ(input_y,z,q);
	  u = B*q;
	}
	
	void computeQ(const VectorXd &input_y,const SXMatrix &z,SXMatrix &q)const{

	  assert_eq(input_y.size(),hatWs.size1());
	  assert_eq(z.size1(),hatWs.size2());
	  assert_eq(z.size2(),1);
	  const SXMatrix in_y = CASADI::convert(input_y);
	  const SXMatrix y = in_y + hatWs.mul(z);
	  computeQ_y(y,q);
	}
	void computeQ(const SXMatrix &z,SXMatrix &q)const{
	  const VectorXd y = VectorXd::Zero(hatWs.size1());
	  computeQ(y,z,q);
	}
	void computeQ(const VectorXd &input_y,const VectorXd &z,VectorXd &q)const{

	  const SXMatrix zs = CASADI::convert(z);
	  const SXMatrix in_y = CASADI::convert(input_y);
	  const SXMatrix y = in_y + hatWs.mul(zs);
	  SXMatrix qs;
	  computeQ_y(y,qs);
	  q = CASADI::convert2Vec<double>(qs);
	}
	void computeQ(const VectorXd &z,VectorXd &q)const{
	  const VectorXd y = VectorXd::Zero(hatWs.size1());
	  computeQ(y,z,q);
	}
	void computeQ_y(const SXMatrix &y,SXMatrix &q)const{
	  SXMatrix b;
	  computeB(y,b);
	  assert_eq(b.size1(),Ps.size2());
	  q = Ps.mul(b);
	}
	void computeQ_y(const VectorXd &y,VectorXd &q)const{
	  const SXMatrix ys = CASADI::convert(y);
	  SXMatrix qs;
	  computeQ_y(ys,qs);
	  q = CASADI::convert2Vec<double>(qs);
	}

	const MatrixXd &getBasisMat()const{
	  return B;
	}
	const MatrixXd &getHatW()const{
	  return hatW;
	}
	const SXMatrix &getHatWs()const{
	  return hatWs;
	}
	const SparseMatrix<double> &get_G()const{
	  return G;
	}
	VectorXd getCubatureY(const VectorXd &full_y)const{

	  VectorXd sub_y(sorted_selPoints.size()*9);
	  for (size_t i = 0; i < sorted_selPoints.size(); ++i){
		assert_in(sorted_selPoints[i]*9,0,(int)full_y.size()-9);
		sub_y.segment(i*9,9) = full_y.segment(sorted_selPoints[i]*9,9);
	  }
	  return sub_y;
	}
	const vector<double> &getSqrtV()const{
	  return Sqrt_V;
	}

  protected:
	void setCubature(const vector<int> &selPoints, const vector<double> &cubWeights){
	  
	  // re-shape matrix
	  SparseMatrix<double> S;
	  set<int> selPointsSet;
	  for (size_t i = 0; i < selPoints.size(); ++i){
	  	assert_in(selPoints[i],0,P.cols()/9);
	  	selPointsSet.insert(selPoints[i]);
	  }
	  EIGEN3EXT::genReshapeMatrix(P.cols(),9,selPointsSet,S,false);
	  
	  // Pj, BPj, hatWj
	  assert_eq(selPoints.size(),cubWeights.size());
	  for (size_t i = 0; i < selPointsSet.size(); ++i){
	  	assert_gt(cubWeights[i],0);
	  	P.block(0,selPoints[i]*9,P.rows(),9) *= cubWeights[i];
	  }
	  P = P*S.transpose();
	  hatW = S*hatW;

	  // Sqrt_Vj
	  sorted_selPoints = selPoints;
	  std::sort(sorted_selPoints.begin(), sorted_selPoints.end());
	  vector<double> tempt(sorted_selPoints.size());
	  for (size_t i = 0; i < sorted_selPoints.size(); ++i){
	  	assert_in(sorted_selPoints[i],0,(int)Sqrt_V.size());
	  	tempt[i] = Sqrt_V[sorted_selPoints[i]];
	  }
	  Sqrt_V = tempt;
	}
	void init(const RS2Euler &rs2euler,const MatrixXd &NLBasis,const MatrixXd &LBasis){

	  this->G = rs2euler.get_G();
	  B = NLBasis;
	  const SparseMatrix<double> &A = rs2euler.get_L();
	  const MatrixXd T = (B.transpose()*(A*B));
	  P = T.inverse();
	  P = P*(B.transpose()*rs2euler.get_VG_t());
	  MapMA2RS::computeMapMatPGW(rs2euler.get_G(),LBasis,hatW);
	  Sqrt_V = rs2euler.get_Sqrt_V();
	}
	void computeB(const SXMatrix &y, SXMatrix &b)const{

	  assert_eq(y.size2(),1);
	  assert_eq(y.size1(),(int)Sqrt_Vs.size()*9);

	  b.resize(y.size1(),1);
	  for (int i = 0; i < y.size1(); i+=9){
		SXMatrix yj(9,1);
		for (int k = 0; k < 9; ++k){
		  yj(k,0) = y(i+k,0);
		}
		yj(0,0) = yj(0,0) + 1e-10;//@todo
		const SX Sqrt_Vj = Sqrt_Vs[i/9];
		const SXMatrix bj = ComputeBjAD::compute(yj,Sqrt_Vj);
		assert_eq(bj.size1(),9);
		assert_eq(bj.size2(),1);
		for (int k = 0; k < 9; ++k){
		   b(i+k,0) = bj(k,0);
		}
	  }
	}
	void convert2Symbols(){

	  CASADI::convert(B,Bs);
	  CASADI::convert(P,Ps);
	  CASADI::convert(hatW,hatWs);
	  Sqrt_Vs.resize(Sqrt_V.size());
	  for (size_t i = 0; i < Sqrt_V.size(); ++i){
		Sqrt_Vs[i] = Sqrt_V[i];
	  }
	}

  private:
	MatrixXd B; // nonlinear basis
	MatrixXd P; // (B^t*A*B)^{-1}B^T*((VG)^t)
	MatrixXd hatW;
	vector<double> Sqrt_V;	

	SXMatrix Bs; // nonlinear basis
	SXMatrix Ps; // (B^t*A*B)^{-1}B^T*((VG)^t)
	SXMatrix hatWs;
	vector<SX> Sqrt_Vs;
	SparseMatrix<double> G;
	vector<int> sorted_selPoints;
  };
  
  /**
   * @class RedRSWarperGrad compute the gradient of each tetrahedron with
   * respect to z using auto-diff. 
   * @warnning: this method is very slow.
   * 
   */
  class RedRSWarperGradAD{
	
  public:
	RedRSWarperGradAD(){
	  
	  vector<CasADi::SX> y = CASADI::makeSymbolic(9,"y");
	  const CasADi::SXMatrix yj = CASADI::convert(y);
	  const CasADi::SX Sqrt_Vj = CasADi::SX("v");
	  const CasADi::SXMatrix bj = ComputeBjAD::compute(yj,Sqrt_Vj);

	  BjFun = CasADi::SXFunction(y,bj);
	  BjFun.init();
	  const CasADi::SXMatrix JacMat = BjFun.jac();
	  y.push_back(Sqrt_Vj);
	  Jacobian = CasADi::SXFunction(y,JacMat);
	  Jacobian.init();
	}
	void computeB(const VectorXd &y,const vector<double>&Sqrt_V,VectorXd &b){

	  assert_eq(y.size(),(int)Sqrt_V.size()*9);
	  b.resize(y.size());
	  const int elemNum = (int)Sqrt_V.size();
	  for (int i = 0; i < elemNum; ++i){
		Eigen::Matrix<double,9,1> bj;
		evaluateFun(Sqrt_V[i],y.segment(i*9,9),bj);
		b.segment(i*9,9) = bj;
	  }
	}
	void jacobian(const VectorXd &y, 
				  const vector<double> &Sqrt_V,
				  const int nodeId, 
				  const MatrixXd &hatW,
				  const MatrixXd &BP,
				  Eigen::Matrix<double,3,-1> &Jz){

	  const int elemNum = y.size()/9;
	  Eigen::Matrix<double,9,9> J;
	  dBdz.resize(y.size(),hatW.cols());
	  for (int i = 0; i < elemNum; ++i){
		evaluateGrad(Sqrt_V[i],y.segment(i*9,9),J);
		dBdz.block(i*9,0,9,hatW.cols()) = J*hatW.block(9*i,0,9,hatW.cols());
	  }
	  Jz = BP.block(nodeId*3,0,3,BP.cols())*dBdz;
	}
	void evaluateGrad(double Sqrt_V, const Eigen::Matrix<double,9,1> &yj, Eigen::Matrix<double,9,9> &rlst){
	  
	  VectorXd y(10);
	  y.segment(0,9) = yj;
	  if (y.segment(0,3).norm() < 1e-8){
		y[0] = 1e-8;
	  }
	  y[9] = Sqrt_V;
	  Jacobian.setInput(&y[0]);
	  Jacobian.evaluate();
	  const CasADi::DMatrix &out = Jacobian.output();
	  assert_eq(out.size1(),9);
	  assert_eq(out.size2(),9);
	  for (int i = 0; i < 9; ++i){
		for (int j = 0; j < 9; ++j){
		  rlst(i,j) = out(i,j).toScalar();
		}
	  }
	}
	void evaluateFun(double Sqrt_V, const Eigen::Matrix<double,9,1> &yj, Eigen::Matrix<double,9,1> &rlst){
	  
	  VectorXd y(10);
	  y.segment(0,9) = yj;
	  if (y.segment(0,3).norm() < 1e-8){
		y[0] = 1e-8;
	  }
	  y[9] = Sqrt_V;
	  BjFun.setInput(&y[0]);
	  BjFun.evaluate();
	  const std::vector<double> &data = BjFun.output().data();
	  assert_eq(data.size(),9);
	  for (int j = 0; j < 9; ++j){
		rlst(j,0) = data[j];
	  }
	}

  private:
	CasADi::SXFunction BjFun; // function for computing Bj.
	CasADi::SXFunction Jacobian; //Jacobian for tetrahedron i with respect to yj.
	MatrixXd dBdz;
  };
  
  typedef boost::shared_ptr<RedRSWarperAD> pRedRSWarperAD;

}//end of namespace


#endif /* _REDRSWARPERAD_H_ */
