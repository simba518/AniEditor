#ifndef _WARPEDPOSCONAD_H_
#define _WARPEDPOSCONAD_H_

#include <boost/shared_ptr.hpp>
#include <RSWarperExt.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <symbolic/sx/sx_tools.hpp>
#include <symbolic/matrix/matrix_tools.hpp>
#include <symbolic/sx/sx.hpp>
#include <symbolic/fx/sx_function.hpp>
#include <symbolic/stl_vector_tools.hpp>
#include <ComputeBj.h>
#include <assertext.h>
using namespace Eigen;
using namespace IEDS;

namespace LSW_ANI_EDITOR{
  
  /**
   * @class WarpedPosConAD compute the matrices for warped position constraints.
   * 
   */
  class WarpedPosConAD{
	
  public:
	WarpedPosConAD(pRSWarperExt warper):_warper(warper){}
	void prepare(const MatrixXd &W);
	void compute(const VectorXd &z_i, const int frame_id, MatrixXd &hat_B, VectorXd &hat_b);
	void get_p2Fp2u(SparseMatrix<double> &p2Fp2u)const;

  protected:
	void compute_B(MatrixXd &B, const VectorXd &r_y);
	void compute_b(VectorXd &b, const VectorXd &r_y);
	const VectorXd compute_r(const VectorXd &z, const VectorXd &y)const;
	CasADi::SXMatrix V9toM3x3(const CasADi::SXMatrix &v)const;
	const VectorXd M3x3toV9(const LSW_WARPING::mat3r &m)const;
	CasADi::SXMatrix assembleS(const CasADi::SXMatrix y_wz)const;
	
  private:
	CasADi::SXFunction _fun_B;
	CasADi::SXFunction _fun_b;
	CasADi::SXFunction _fun_pFpu;
	CasADi::SXMatrix _p2Fp2u;

	SparseMatrix<double> _inv_A1;
	SparseMatrix<double> _inv_A2;
	pRSWarperExt _warper;
	MatrixXd eigen_hat_W; // warp modal matrix
  };
  
  typedef boost::shared_ptr<WarpedPosConAD> pWarpedPosConAD;
  
}//end of namespace

#endif /*_WARPEDPOSCONAD_H_*/
