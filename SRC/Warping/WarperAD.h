#ifndef _WARPERAD_H_
#define _WARPERAD_H_

#include <string.h>
#include <CASADITools.h>
#include <Timer.h>
#include <DefGradOperator.h>
#include <TetMesh.h>

using namespace CASADI;
using UTILITY::Timer;
using CasADi::SXMatrix;
using CasADi::SX;
using CasADi::SXFunction;
using namespace UTILITY;

namespace LSW_WARPING{
  
  /**
   * @class WarperAD given linear displacement p, compute warped displacement u.
   * 
   */
  class WarperAD{
	
  public:
	static void computeWarpFun(pTetMesh_const rest, 
							   const SXMatrix &p,
							   SXMatrix &u){
	  vector<SXMatrix> w;
	  computeRotVector(rest,p,w);
	  computeWarpFun(p,w,u);
	}

	static void computeWarpFun(const SXMatrix &p,
							   const vector<SXMatrix> &nodeW,
							   SXMatrix &u){
	  
	  u.resize(p.size1(),1);
	  SXMatrix v3(3,1);
	  for (size_t i = 0; i < nodeW.size(); ++i){
		const SXMatrix warpR = warpMat(nodeW[i]);
		v3(0) = p(i*3+0);
		v3(1) = p(i*3+1);
		v3(2) = p(i*3+2);
		v3 = warpR.mul(v3);
		u(i*3+0) = v3(0);
		u(i*3+1) = v3(1);
		u(i*3+2) = v3(2);
	  }
	}

	static void computeRotMat(pTetMesh_const rest, const MatrixXd &W, MatrixXd &rotW){
	  
	  assert(rest);
	  assert_eq(rest->nodes().size()*3,W.rows());
	  SparseMatrix<double> rotM;
	  computeRotMat(rest,rotM);
	  assert_eq(rest->nodes().size()*3,rotM.cols());
	  assert_eq(rest->nodes().size()*3,rotM.rows());
	  rotW = rotM*W;
	}

	static void computeRotMat(pTetMesh_const rest, SparseMatrix<double> &rotM){
	  
	  const vector<SX> p = CASADI::makeSymbolic(rest->nodes().size()*3,"p");
	  const SXMatrix ps = CASADI::convert(p);
	  vector<SXMatrix> nodeW;
	  computeRotVector(rest,ps,nodeW);
	  SXMatrix w(nodeW.size()*3,1);
	  for (size_t i = 0; i < nodeW.size(); ++i){
		w(i*3+0) = nodeW[i](0);
		w(i*3+1) = nodeW[i](1);
		w(i*3+2) = nodeW[i](2);
	  }
	  SXFunction wp = SXFunction(p,w);
	  wp.init();
	  const SXMatrix rotMs = wp.jac();
	  CASADI::convert(rotMs,rotM);
	}

	static void computeRotVector(pTetMesh_const rest, 
								 const SXMatrix &p,
								 vector<SXMatrix> &nodeW){
	  assert(rest);
	  assert_eq(rest->nodes().size(),p.size1()/3);
	  assert_eq(p.size1()%3,0);
	  assert_eq(p.size2(),1);

	  // compute gradient operator G
	  SparseMatrix<double> G;
	  DefGradOperator::compute(rest,G);
	  SXMatrix sG;
	  CASADI::convert(G,sG);
	  const SXMatrix m = mul(sG,p);
	  const int elem_num = m.size()/9;
	  assert_eq(rest->tets().size(),elem_num);

	  // compute w for each element
	  vector<SXMatrix> elemW(elem_num);
	  SXMatrix sw(3,1);
	  for (int j =0; j < elem_num; ++j){
	  	sw(0) = (m(j*9+7)-m(j*9+5))/2;
	  	sw(1) = (m(j*9+2)-m(j*9+6))/2;
	  	sw(2) = (m(j*9+3)-m(j*9+1))/2;
	  	elemW[j] = sw;
	  }

	  // compute w for each node
	  vector<int> eleCount(p.size1()/3);
	  nodeW.resize(p.size1()/3);

	  sw.setZero();
	  for (size_t i = 0; i < nodeW.size(); ++i){
		nodeW[i] = sw;
		eleCount[i] = 0;
	  }
	  for (int j = 0; j < elem_num; ++j){
		for (int v = 0; v < 4; ++v){
		  const int i = rest->tets()[j][v];
		  assert_in(i,0,(int)nodeW.size()-1);
		  nodeW[i] += elemW[j];
		  eleCount[i] ++;
		}
	  }
	  for (size_t i = 0; i < nodeW.size(); ++i){
		assert_gt (eleCount[i],0);
		nodeW[i] /= (double)eleCount[i];
	  }
	}

	static SXMatrix warpMat(const SXMatrix &w){

	  assert_eq(w.size1(),3);
	  assert_eq(w.size2(),1);
	  const SXMatrix nw = norm_2(w);
	  const SXMatrix I  = SXMatrix::eye(3);
	  const SXMatrix _w = w/nw;
	  const SXMatrix wx = SkewMat(_w);
	  return I + wx*( (1.0f-cos(nw))/nw ) + (mul(wx,wx))*( 1.0f-sin(nw)/nw );
	}

	static SXMatrix SkewMat(const SXMatrix &_w){

	  assert_eq(_w.size1(),3);
	  assert_eq(_w.size2(),1);

	  SXMatrix wx(3,3);
	  wx.setZero();
	  
	  wx(0,1) = -_w[2];
	  wx(0,2) = _w[1];
	  wx(1,2) = -_w[0];

	  wx(1,0) = -wx(0,1);
	  wx(2,0) = -wx(0,2);
	  wx(2,1) = -wx(1,2);

	  return wx;
	}
  };

  /**
   * @class WarperExtAD compute the deformation Jacobian of u(p+Wz) with respect
   * to z, where p is the input local (linear) displacement, and W is the modal
   * matrix, while z is the modal coordinates.
   * 
   */
  class WarperExtAD{
	
  public:
	void precompute(pTetMesh_const rest, const MatrixXd &W){

	  const vector<SX> inputP = CASADI::makeSymbolic(W.rows(),"p");
	  const vector<SX> z = CASADI::makeSymbolic(W.cols(),"z");
	  SXMatrix u;
	  computeWarpFun(rest,W,z,inputP,u);
	  vector<SX> x = inputP;
	  x.insert(x.end(),z.begin(),z.end());
	  SXFunction uz = SXFunction(z,u);
	  uz.init();
	  const SXMatrix JacMat = uz.jac();
	  Jacobian = SXFunction(x,JacMat);
	  Jacobian.init();
	}
	void jacobian(const VectorXd &inputP,const VectorXd &z, MatrixXd &J)const{

	  VectorXd x(inputP.size()+z.size());
	  x.segment(0,inputP.size()) = inputP;
	  x.segment(inputP.size(),z.size()) = z;
	  SXFunction jac = Jacobian;
	  CASADI::evaluate(jac,x,J);
	}

  protected:
	void computeWarpFun(pTetMesh_const rest, const MatrixXd &W, 
						const vector<SX> &z,
						const vector<SX> &inputP,
						SXMatrix &u){
	  
	  assert(rest);
	  assert_gt(W.size(),0);
	  assert_eq(rest->nodes().size()*3,W.rows());

	  SXMatrix zs(z.size(),1), ps(inputP.size(),1);
	  for (size_t i = 0; i < z.size(); ++i)
		zs(i) = z[i];
	  for (size_t i = 0; i < inputP.size(); ++i)
		ps(i) = inputP[i];

	  const SXMatrix sW = CASADI::convert(W);
	  const SXMatrix p = ps + sW.mul(zs);
	  WarperAD::computeWarpFun(rest,p,u);
	}
	
  private:
	SXFunction Jacobian;
  };

  /**
   * @class WarperPerNodeAD compute the deformation Jacobian of each node i,
   * i.e., u_i(w_i,p_i,z) with respect to z, where p_i is the input local
   * (linear) displacement of node i, w_i is the input rotation vector of node
   * i, while z is the modal coordinates.
   * 
   */
  class WarperPerNodeAD{
	
  public:
	void precompute(pTetMesh_const rest, const MatrixXd &W){

	  // prepare symbols
	  const vector<SX> inputW = CASADI::makeSymbolic(W.rows(),"w");
	  const vector<SX> inputP = CASADI::makeSymbolic(W.rows(),"p");
	  const vector<SX> z = CASADI::makeSymbolic(W.cols(),"z");

	  // compute p = inputP+W*z
	  const SXMatrix zs = CASADI::convert(z);
	  const SXMatrix ps = CASADI::convert(inputP);
	  const SXMatrix Ws = CASADI::convert(W);
	  const SXMatrix Wz = Ws.mul(zs);
	  const SXMatrix p = ps + Wz;

	  // compute w = inputW + w(z)
	  vector<SXMatrix> nodeW;
	  WarperAD::computeRotVector(rest,Wz,nodeW);
	  for (size_t i = 0; i < nodeW.size(); ++i){
		nodeW[i](0) += inputW[i*3+0];
		nodeW[i](1) += inputW[i*3+1];
		nodeW[i](2) += inputW[i*3+2];
	  }

	  // compute u
	  SXMatrix u;
	  WarperAD::computeWarpFun(p,nodeW,u);

	  // compute grad(u[i]) for each node i
	  Jacobian.resize(nodeW.size());
	  SXMatrix u3(3,1);
	  for (size_t i = 0; i < Jacobian.size(); ++i){

		// assemble u_i(z)
		u3(0) = u(i*3+0);
		u3(1) = u(i*3+1);
		u3(2) = u(i*3+2);
		SXFunction uz = SXFunction(z,u3);
		uz.init();

		// assemble function for d(u_i)/d(z)(x), for x = [w,p,z].
		vector<SX> x;
		x.insert(x.end(),inputW.begin()+i*3,inputW.begin()+i*3+3);
		x.insert(x.end(),inputP.begin()+i*3,inputP.begin()+i*3+3);
		x.insert(x.end(),z.begin(),z.end());
		const SXMatrix JacMat = uz.jac();
		Jacobian[i] = SXFunction(x,JacMat);
		Jacobian[i].init();
	  }
	}
	void jacobian(const Vector3d &inputW3, const Vector3d &inputP3, 
				  const VectorXd &z, const int nodeId, MatrixXd &J)const{

	  assert_in(nodeId,0,Jacobian.size()-1);
	  VectorXd x(z.size()+6); // [w3,p3,z]
	  x.segment(0,3) = inputW3;
	  x.segment(3,3) = inputP3;
	  x.segment(6,z.size()) = z;
	  SXFunction jac = Jacobian[nodeId];
	  CASADI::evaluate(jac,x,J);
	}
	
  private:
	vector<SXFunction> Jacobian; // Jacobian[i] is the Jacobian for node i.
  };

  /**
   * @class WarperPerNodeAD_WP compute the deformation Jacobian of each node i,
   * i.e., u_i(w_i,p_i) with respect to w_i and p_i respectively, where p_i is
   * the local (linear) displacement of node i, w_i is the rotation vector of
   * node i.
   * 
   */
  class WarperPerNodeAD_WP{
	
  public:
	void precompute(pTetMesh_const rest){
  
	  // prepare symbols
	  const vector<SX> w = CASADI::makeSymbolic(rest->nodes().size()*3,"w");
	  const vector<SX> p = CASADI::makeSymbolic(rest->nodes().size()*3,"p");
	  const SXMatrix ps = CASADI::convert(p);
	  const vector<SXMatrix> nodeW = CASADI::Sx2Mat(w,w.size()/3);

	  // compute u
	  SXMatrix u;
	  WarperAD::computeWarpFun(ps,nodeW,u);

	  // compute d(u[i])/dw[i] and  d(u[i])/dp[i] for each node i
	  JacobianW.resize(nodeW.size());
	  JacobianP.resize(nodeW.size());
	  SXMatrix u3(3,1);
	  vector<SX> w3(3),p3(3);
	  for (size_t i = 0; i < JacobianP.size(); ++i){

		// assemble u_i(z) and x = [wi,pi]
		u3(0) = u(i*3+0);
		u3(1) = u(i*3+1);
		u3(2) = u(i*3+2);
		vector<SX> x(6);
		x[0] = w[i*3+0];		x[1] = w[i*3+1];		x[2] = w[i*3+2];
		x[3] = p[i*3+0];		x[4] = p[i*3+1];		x[5] = p[i*3+2];

		// assemble u(w), compute du/dw
		w3[0] = w[i*3+0]; 
		w3[1] = w[i*3+1]; 
		w3[2] = w[i*3+2]; 
		SXFunction uw = SXFunction(w3,u3);
		uw.init();
		const SXMatrix JacMatW = uw.jac();
		JacobianW[i] = SXFunction(x,JacMatW);
		JacobianW[i].init();

		// assemble u(p), compute du/dw
		p3[0] = p[i*3+0]; 
		p3[1] = p[i*3+1]; 
		p3[2] = p[i*3+2]; 
		SXFunction up = SXFunction(p3,u3);
		up.init();
		const SXMatrix JacMatP = up.jac();
		JacobianP[i] = SXFunction(x,JacMatP);
		JacobianP[i].init();
	  }
	}
	void jacobian(const Vector3d &w3, const Vector3d &p3, const int nodeId, 
				   Matrix3d &Jw, Matrix3d &Jp){

	  assert_in(nodeId,0,JacobianW.size()-1);
	  assert_in(nodeId,0,JacobianP.size()-1);
	  Matrix<double, 6, 1 > x;
	  (x.block<3,1>)(0,0) = w3;
	  (x.block<3,1>)(3,0) = p3;
	  evaluate(JacobianW[nodeId],x,Jw);
	  evaluate(JacobianP[nodeId],x,Jp);
	}

  protected:
	void evaluate(SXFunction &fun, const Matrix<double,6,1> &x, Matrix3d &rlst)const{
	  
	  fun.setInput(&x[0]);
	  fun.evaluate();
	  const std::vector<double> &data = fun.output().data();
	  assert_eq(data.size(),9);
	  rlst(0,0) = data[0];
	  rlst(0,1) = data[1];
	  rlst(0,2) = data[2];

	  rlst(1,0) = data[3];
	  rlst(1,1) = data[4];
	  rlst(1,2) = data[5];

	  rlst(2,0) = data[6];
	  rlst(2,1) = data[7];
	  rlst(2,2) = data[8];
	}
	
  private:
	vector<SXFunction> JacobianW; // Jacobian[i] is the Jacobian for node i with respect to wi.
	vector<SXFunction> JacobianP; // Jacobian[i] is the Jacobian for node i with respect to pi.
  };

}//end of namespace

#endif /*_WARPERAD_H_*/
