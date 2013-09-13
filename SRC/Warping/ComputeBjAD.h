#ifndef _COMPUTEBJAD_H_
#define _COMPUTEBJAD_H_

#include <vector>
#include <symbolic/sx/sx.hpp>
#include <symbolic/matrix/matrix_tools.hpp>
#include <assertext.h>

namespace LSW_WARPING{
  
  /**
   * @class ComputeBjAD compute Bj using AD.
   * @see ComputeBj
   */
  class ComputeBjAD{
	
  public:
	// compute bj
	static void compute(const CasADi::SXMatrix &yj,CasADi::SX Sqrt_Vj, CasADi::SXMatrix &bj){
	  
	  assert_eq(yj.size1(),9);
	  assert_eq(yj.size2(),1);
	  
	  // the symtric part
	  CasADi::SXMatrix SI(3,3);
	  compute_S(yj, SI);
	  SI(0,0) += 1.0f;
	  SI(1,1) += 1.0f;
	  SI(2,2) += 1.0f;
  
	  // the antisymmetric part
	  CasADi::SXMatrix Exp_W(3,3);
	  compute_exp(yj, Exp_W);

	  // bj
	  CasADi::SXMatrix m = Sqrt_Vj*(Exp_W.mul(SI));
	  m(0,0) -= Sqrt_Vj;
	  m(1,1) -= Sqrt_Vj;
	  m(2,2) -= Sqrt_Vj;

	  bj.resize(9,1);
	  bj(0) = m(0,0);  bj(1) = m(0,1);  bj(2) = m(0,2);// row 0
	  bj(3) = m(1,0);  bj(4) = m(1,1);  bj(5) = m(1,2);// row 1
	  bj(6) = m(2,0);  bj(7) = m(2,1);  bj(8) = m(2,2);// row 2
	}

	static CasADi::SXMatrix compute(const CasADi::SXMatrix &yj,CasADi::SX Sqrt_Vj){
	  CasADi::SXMatrix b;
	  compute(yj,Sqrt_Vj,b);
	  return b;
	}

  protected:
	static const CasADi::SXMatrix &compute_S(const CasADi::SXMatrix &yj, CasADi::SXMatrix &SI){

	  SI.resize(3,3);
	  SI(0,0) = yj(0+3);
	  SI(1,0) = yj(1+3);
	  SI(2,0) = yj(2+3);
	  SI(1,1) = yj(3+3);
	  SI(2,1) = yj(4+3);
	  SI(2,2) = yj(5+3);

	  SI(0,1) = SI(1,0);
	  SI(0,2) = SI(2,0);
	  SI(1,2) = SI(2,1);
	  return SI;
	}

	static const CasADi::SXMatrix &compute_exp(const CasADi::SXMatrix &yj, CasADi::SXMatrix &Exp_W){

	  Exp_W.resize(3,3);
	  const std::vector<CasADi::SX> &dy = yj.data();
	  assert_eq(dy.size(),9);
	  const CasADi::SX theta = sqrt(dy[0]*dy[0] + dy[1]*dy[1] + dy[2]*dy[2]);
	  const CasADi::SX x = dy[0]/theta;
	  const CasADi::SX y = dy[1]/theta;
	  const CasADi::SX z = dy[2]/theta;
	  const CasADi::SX c = 1.0f - cos(theta);
	  Exp_W(0,0) = cos(theta)+x*x*c ;      Exp_W(0,1) = x*y*c-z*sin(theta);     Exp_W(0,2) = x*z*c + y* sin(theta);  // row 0
	  Exp_W(1,0) = x*y *c + z*sin(theta) ; Exp_W(1,1) = cos(theta) + y*y *c ;   Exp_W(1,2) =  y* z*c - x*sin(theta); // row 1
	  Exp_W(2,0) = z*x*c - y* sin(theta) ; Exp_W(2,1) = z*y *c + x*sin(theta) ; Exp_W(2,2) =  cos(theta) + z*z *c;   // row 2
	  return Exp_W;
	}

  };
  
}//end of namespace

#endif /*_COMPUTEBJAD_H_*/
