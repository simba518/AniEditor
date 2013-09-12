#include <boost/lexical_cast.hpp>
#include <symbolic/sx/sx.hpp>
#include <symbolic/sx/sx_tools.hpp>
#include "ComputeGeAutoDiff.h"
#include <assertext.h>
using namespace LSW_WARPING;

pComputeGeAutoDiff ComputeGeAutoDiff::p_instance;

void ComputeGeAutoDiff::initFun(){
  
  CasADi::SXMatrix P(3,4), new_P(3,4);
  for (int i = 0; i < 3; ++i){
	for (int j = 0; j < 4; ++j){
	  const string si = boost::lexical_cast<string>(i);
	  const string sj = boost::lexical_cast<string>(j);
	  P.elem(i,j) = CasADi::SX("v_"+si+sj);
	  new_P.elem(i,j) = CasADi::SX("new_v"+si+sj);
	} 
  }

  CasADi::SXMatrix V(3,3), new_V(3,3);
  for (int c = 0; c < 3; ++c){
	for (int r = 0;  r < 3; ++r){
	  V.elem(r,c) = P.elem(r,c+1)-P.elem(r,0);
	  new_V.elem(r,c) = new_P.elem(r,c+1)-new_P.elem(r,0);
	}
  }

  const CasADi::SXMatrix M = new_V.mul(CasADi::inv(V));

  vector<CasADi::SX> v(12), new_v(12), all_v(24);
  for (int i = 0; i < 3; ++i){
	for (int j = 0; j < 4; ++j){
	  v[j*3+i] = P.elem(i,j);
	  new_v[j*3+i] = new_P.elem(i,j);

	  all_v[j*3+i] = P.elem(i,j);
	  all_v[12+j*3+i] = new_P.elem(i,j);
	} 
  }

  CasADi::SXMatrix m = CasADi::SXMatrix::zeros(9,1);
  m[0] = M.elem(0,0);   m[1] = M.elem(0,1);   m[2] = M.elem(0,2); // row 0
  m[3] = M.elem(1,0);   m[4] = M.elem(1,1);   m[5] = M.elem(1,2); // row 1
  m[6] = M.elem(2,0);   m[7] = M.elem(2,1);   m[8] = M.elem(2,2); // row 2

  const CasADi::SXMatrix Ge = CasADi::jacobian(m,new_v);
  Ge_fun = CasADi::SXFunction(v,Ge);
  Ge_fun.init();
  Ge.sparsity().getSparsity (rows, cols);
  assert_eq (cols.size(),rows.size());
}
