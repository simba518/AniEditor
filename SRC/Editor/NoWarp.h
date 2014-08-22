#ifndef _NOWARP_H_
#define _NOWARP_H_

#include <boost/shared_ptr.hpp>

namespace LSW_ANI_EDITOR{
  
  /**
   * @class NoWarp we construct the displacements using u(z) = W*z without
   * warping, and here we give the function to compute the function and
   * gradient.
   * 
   */
  class NoWarp{
	
  public:
	NoWarp(const MatrixXd &W){
	  this->W = W;
	}
	void warp(const VectorXd &z,int frame_id,const vector<int> &nodes,VectorXd &ui){
	  assert_eq(z.size(),W.cols());
	  ui = block(W,nodes)*z;
	}
	void jacobian(const VectorXd &z,int frame_id,const vector<int> &nodes,
				  const VectorXd &uc,VectorXd &g){
	  const MatrixXd Wb = block(W,nodes);
	  g = Wb.transpose()*(Wb*z-uc);
	}
	CasADi::SXMatrix warp(const CasADi::SXMatrix&z,const vector<int> &nodes){
	  assert_eq(z.size1(),W.cols());
	  assert_eq(z.size2(),1);
	  CasADi::SXMatrix sWi;  
	  CASADI::convert(this->block(W,nodes),sWi);
	  return sWi.mul(z);
	}
  
  protected:
	MatrixXd block(const MatrixXd &W,const vector<int> &nodes)const{

	  MatrixXd Wi(nodes.size()*3, W.cols());
	  for (size_t i = 0; i < nodes.size(); ++i){
		const int j3 = nodes[i]*3;
		assert_le(j3,W.rows()-3);
		Wi.block( i*3,0,3, W.cols() ) = W.block( j3,0,3,W.cols() );
	  }
	  return Wi;
	}
  
  private:
	MatrixXd W;
  };
  
  typedef boost::shared_ptr<NoWarp> pNoWarp;
  
}//end of namespace

#endif /*_NOWARP_H_*/
