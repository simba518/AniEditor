#include <symbolic/matrix/slice.hpp>
#include <CASADITools.h>
#include <SparseMatrixTools.h>
#include <MapMA2RS.h>
#include "WarpedPosConAD.h"
using namespace LSW_ANI_EDITOR;

void WarpedPosConAD::prepare(const MatrixXd &W){

  // prepare A, V, G, 
  const SparseMatrix<double> &A = _warper->getRS2Euler().get_A();
  const vector<double> &Sqrt_V = _warper->getRS2Euler().get_Sqrt_V();
  const SparseMatrix<double> &eigen_G = _warper->get_G();
  CasADi::SXMatrix G;
  CASADI::convert(eigen_G,G);

  // compute wapred(W): hat_W.
  LSW_WARPING::MapMA2RS::computeMapMatPGW(eigen_G, W, eigen_hat_W);
  CasADi::SXMatrix hat_W;
  CASADI::convert(eigen_hat_W, hat_W);

  // get dimensions
  const int full_dim = W.rows();
  const int red_dim = W.cols();
  const int tet_numx9 = eigen_hat_W.rows();
  const int tet_num = tet_numx9/9;

  // u, r, y, z
  const CasADi::SXMatrix u = CasADi::ssym("u",full_dim);
  const CasADi::SXMatrix z = CasADi::ssym("z",red_dim);
  const CasADi::SXMatrix r = CasADi::ssym("r",tet_numx9);
  const CasADi::SXMatrix y = CasADi::ssym("y",tet_numx9);
  CasADi::SXMatrix u_z = u;   u_z.append(z);
  CasADi::SXMatrix r_y = r;   r_y.append(y);
  CasADi::SXMatrix x = u_z;   x.append(r_y);

  // assemble objective function: F
  CasADi::SXMatrix F = 0.0f;
  const CasADi::SXMatrix wz = hat_W.mul(z);
  const CasADi::SXMatrix gu = G.mul(u);
  const CasADi::SXMatrix I = CasADi::SXMatrix::eye(3);
  for (int j = 0; j < tet_num; ++j){

  	const CasADi::Slice seg9(j*9,j*9+9,1);
  	const CasADi::SXMatrix R  = V9toM3x3( r.getSub(seg9,0));
  	const CasADi::SXMatrix Gu = V9toM3x3(gu.getSub(seg9,0));
  	const CasADi::SXMatrix S  = assembleS( y.getSub(seg9,0) + wz.getSub(seg9,0) );
  	const CasADi::SXMatrix F1 = (I + Gu - R.mul(I+S))*Sqrt_V[j];
  	F += 0.5*CasADi::sumAll(F1*F1);
  }

  // make functions.
  CasADi::SXMatrix pFpu = CasADi::jacobian(F,u);
  _fun_pFpu = CasADi::SXFunction(x,pFpu);
  _fun_pFpu.init();

  CasADi::SXMatrix u0z0_ry = CasADi::SXMatrix::zeros(u_z.size(),1);
  u0z0_ry.append(r_y);
  const CasADi::SXMatrix b = _fun_pFpu.eval(u0z0_ry);

  _fun_b = CasADi::SXFunction(r_y,b);
  _fun_b.init();

  const CasADi::SXMatrix p2Fpupz = CasADi::jacobian(pFpu,z);
  _fun_B = CasADi::SXFunction(r_y,p2Fpupz);
  _fun_B.init();

  _p2Fp2u = CasADi::jacobian(pFpu,u);

  // compute inv_A1, inv_A2
  SparseMatrix<double> inv_A;
  EIGEN3EXT::inverse(A,inv_A);
  _inv_A1 = EIGEN3EXT::block(inv_A,0,0,u.size(),u.size());
  _inv_A2 = EIGEN3EXT::block(inv_A,0,u.size(),u.size(),A.cols()-u.size());
}

void WarpedPosConAD::compute(const VectorXd &z_i, const int frame_id, MatrixXd &hat_B, VectorXd &hat_b){

  const VectorXd &y = _warper->getInputRSCoord(frame_id);
  const VectorXd r = compute_r(z_i, y);
  VectorXd r_y (r.size() + y.size());
  r_y.segment(0,r.size()) = r;
  r_y.segment(r.size(),y.size()) = y;

  MatrixXd B;
  VectorXd b;
  compute_B(B, r_y);
  compute_b(b, r_y);

  const VectorXd s = VectorXd::Zero(_inv_A2.cols()); /// @bug, use fixed constraints only.
  hat_B = (-1.0)*_inv_A1*B;
  hat_b = (-1.0)*_inv_A1*b + _inv_A2*s;
}

const VectorXd WarpedPosConAD::compute_r(const VectorXd &z, const VectorXd &y)const{
  
  const int tet_num = y.size()/9;
  VectorXd r(y.size());
  LSW_WARPING::mat3r exp_m;
  const VectorXd wz = eigen_hat_W*z;

  for (int i = 0; i < tet_num; ++i) {

	const VectorXd yi = y.segment(9*i,9) + wz.segment(9*i,9);
	LSW_WARPING::ComputeBj::compute_exp(yi, exp_m);
	r.segment(9*i,9) = M3x3toV9(exp_m);
  }
  return r;
}

void WarpedPosConAD::compute_B(MatrixXd &B, const VectorXd &r_y) {

  CASADI::evaluate (_fun_B, r_y, B);
}

void WarpedPosConAD::compute_b(VectorXd &b, const VectorXd &r_y){

  CASADI::evaluate (_fun_b, r_y, b);
}

CasADi::SXMatrix WarpedPosConAD::assembleS(const CasADi::SXMatrix y_wz)const{

  assert_eq (y_wz.size(),9);

  CasADi::SXMatrix S = CasADi::SXMatrix::zeros(3,3);
  S(0,0) = y_wz[0+3];
  S(1,0) = y_wz[1+3];
  S(2,0) = y_wz[2+3];
  S(1,1) = y_wz[3+3];
  S(2,1) = y_wz[4+3];
  S(2,2) = y_wz[5+3];

  S(0,1) = S(1,0);
  S(0,2) = S(2,0);
  S(1,2) = S(2,1);
  return S;
}

void WarpedPosConAD::get_p2Fp2u(SparseMatrix<double> &p2Fp2u)const{
  
  CASADI::convert(_p2Fp2u,p2Fp2u);
}

CasADi::SXMatrix WarpedPosConAD::V9toM3x3(const CasADi::SXMatrix &v)const{

  assert_eq (v.size(),9);
  CasADi::SXMatrix M = CasADi::SXMatrix::zeros(3,3);
  M(0,0) = v(0);      M(0,1) = v(1);      M(0,2) = v(2); // row 0
  M(1,0) = v(3);      M(1,1) = v(4);      M(1,2) = v(5); // row 1
  M(2,0) = v(6);      M(2,1) = v(7);      M(2,2) = v(8); // row 2
  return M;
}

const VectorXd WarpedPosConAD::M3x3toV9(const LSW_WARPING::mat3r &m)const{

  VectorXd v(9);
  memcpy (&v[0], &m(0,0), 9*sizeof(double));
  return v;
}
