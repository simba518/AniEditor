#include <Log.h>
#include "SpaceTimeHessianAutoDiff.h"
using namespace IEDS;

SpaceTimeHessianAutoDiff::SpaceTimeHessianAutoDiff():BaseSpaceTimeHessian(){

}

bool SpaceTimeHessianAutoDiff::compute_HTriplet(vector<E_Triplet> &H_triplet){
 
  makeSymbol_H();
  calculate_H(H_triplet);
  return true;
}

void SpaceTimeHessianAutoDiff::makeSymbol_H(){
  
  TRACE_FUN();
  
  const int r = reducedDim();
  const int T = totalFrames();
  const int n = totalDim();
  assert_ge(T, 3);
  assert_ge(h,0.0f);

  // initialize the symbols of the unknows: z[i]
  vector<CasADi::SX> z(n);
  make_symbolic(z.begin(), z.end(), "z");
  
  // compute the energy function
  const double _h = 1.0f/h;
  E = 0.0f;
  for (int k = 0; k < r; ++k){ 
	
	// compute the energy function of mode k
	CasADi::SX Ek = 0.0f;
    for (int i = 1; i < T-1; ++i){

	  CasADi::SX y0 = z[(i-1)*r + k];
	  CasADi::SX y1 = z[(i+0)*r + k];
	  CasADi::SX y2 = z[(i+1)*r + k];
	  CasADi::SX y = y2 - 2*y1 + y0;
	  const double la = eigen_values[k];
	  CasADi::SX s = (_h*_h*y + _h*(alpha_k*la + alpha_m)*(y2-y1) + la*y1);
	  // CasADi::SX s = (_h*_h*y + _h*(alpha_k*la + alpha_m)*(y2-y0)*0.5f + la*y1);
	  Ek += 0.5*s*s;
	}
	E += Ek;
  }

  // init the function.
  E_z = CasADi::SXFunction(z,E);
  E_z.init();

  // compute the second gradient of E, e.g symbol_H
  H_SX = E_z.hess();
}

void SpaceTimeHessianAutoDiff::calculate_H(vector<E_Triplet> &H_triplet){

  const int n = totalDim();
  assert_ge(H_SX.size1(),n );
  assert_ge(H_SX.size2(),n );

  vector<int> rows,cols;
  H_SX.sparsity().getSparsity (rows, cols);
  const std::vector<CasADi::SX> &data = H_SX.data ();
  
  const int nz = (int)data.size();
  H_triplet.clear();
  H_triplet.reserve(nz);
  for (int i = 0; i < nz; ++i){

	assert_in (rows[i],0,n-1);
	assert_in (cols[i],0,n-1);
	if (data[i].getValue() != 0.0f){
	  H_triplet.push_back( E_Triplet(rows[i],cols[i],data[i].getValue()) );
	}
  }
}
