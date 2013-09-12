#ifndef _COMPUTEGEAUTODIFF_H_
#define _COMPUTEGEAUTODIFF_H_

#include <vector>
#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Sparse>
#include <symbolic/fx/sx_function.hpp>
#include <assertext.h>
using namespace std;

namespace LSW_WARPING{
  
  /**
   * @class ComputeGeAutoDiff compute the discrete deformation gradient operator
   * Ge of one tetrahedron element using auto-differential technology.
   * 
   * @note the gradient is computing using m = Ge*V, where m is a 9x1 vector,
   * which expressed as a row-major 3x3 matrix, that is 
   *                                                        |m0, m1, m2|
   * m = (m0,m1,m2,m3,m4,m5,m6,m7,m8), then the matrix is   |m3, m4, m5|
   *                                                        |m6, m7, m8|
   */
  class ComputeGeAutoDiff{

	typedef Eigen::Triplet<double> T;
	
  public:
	static inline boost::shared_ptr<ComputeGeAutoDiff> getInstance(){
	  
	  if(p_instance == NULL){
		p_instance=boost::shared_ptr<ComputeGeAutoDiff>(new ComputeGeAutoDiff());
	  }
	  return p_instance;
	}

	/**
	 * compute Ge using the vertices in V, where V[i] is the i-th vertice of the
	 * element.
	 * 
	 * @note in the first time before calling this function, Ge[9][12] should be
	 * initialized to zeros.
	 */
	void compute(const double V[4][3], double Ge[9][12]){

	  // evaluate
	  Ge_fun.setInput(&(V[0][0]));
	  Ge_fun.evaluate();
	  const vector<double> &data = Ge_fun.output().data();
	  const int nz = (int)data.size();
	  assert_eq (data.size(),rows.size());
	  
	  // set the nonzeros
	  for (int i = 0; i < nz; ++i){
		const int r = rows[i];
		const int c = cols[i];
		Ge[r][c] = data[i];
	  }
	}
	
  protected:
    ComputeGeAutoDiff(){
	  initFun();
	}

	/// initialize Ge_fun
	void initFun();
	
  private:
	static boost::shared_ptr<ComputeGeAutoDiff> p_instance;
	CasADi::SXFunction Ge_fun; // matrix function for calculating Ge

	// the sparse structure of Ge, where (rows[i],cols[i]) is the nonzero
	// i-th element.
	vector<int> rows; 
	vector<int> cols; 
  };
  
  typedef boost::shared_ptr<ComputeGeAutoDiff> pComputeGeAutoDiff;
  
}//end of namespace

#endif /*_COMPUTEGEAUTODIFF_H_*/
