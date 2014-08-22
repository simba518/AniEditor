#ifndef _BASISSELREGENERGY_H_
#define _BASISSELREGENERGY_H_

#include "MtlOptDM.h"

namespace MTLOPT{
  
  /**
   * @class BasisSelRegEnergy regularization item for the basis selection
   * matrix S.
   * 
   */
  class BasisSelRegEnergy{
	
  public:
	BasisSelRegEnergy(){
	  penalty = 1.0f;
	}

	void setPenalty(const double p){
	  assert_ge(p,0.0f);
	  WARN_LOG_COND("the penalty items for PosConEnergyS is zero.",p!=0);
	  penalty  = p;
	}

	virtual double fun(const double *pS, const int rows, const int cols)const{

	  assert(pS);
	  assert_ge(rows, cols);

	  const MatrixXd &S = Map<MatrixXd>( const_cast<double*>(pS), rows, cols);
	  double f = 0.0f;
	  for (int i = 0; i < S.rows(); ++i){
		const double n = S.row(i).norm();
		f += n*n;
	  }
	  return f*0.5f*penalty;
	}

	virtual void grad(const double *pS, double *pG, const int rows, const int cols)const{

	  assert(pS);
	  assert(pG);
	  assert_ge(rows, cols);

	  const MatrixXd &S = Map<MatrixXd>( const_cast<double*>(pS), rows, cols );

	  static MatrixXd g;
	  g.resize(S.rows(),S.cols());
	  for (int i = 0; i < S.rows(); ++i){
	  	g.row(i) = S.row(i);
	  }
	  g *= penalty;
	  
	  memcpy(pG, &g(0,0), g.size()*sizeof(double) );
	}

	virtual void hessian(VectorXd &H,const int rw, const int rs)const{
	  // hessisan is constant and diagonal.

	  assert_gt(rs,0);
	  assert_ge(rw,rs);

	  H.resize(rs*rw);
	  for (int i = 0; i < rs; ++i)
		H.segment(i*rw,rw).setConstant(1.0f);
	  H *= penalty;
	}

	double getPenalty()const{
	  return penalty;
	}

  protected:
	double penalty;
  };
  typedef boost::shared_ptr<BasisSelRegEnergy> pBasisSelRegEnergy;

  class BasisNormalizeEnergy:public BasisSelRegEnergy{
	
  public:
	double fun(const double *pS, const int rows, const int cols)const{

	  assert(pS);
	  assert_ge(rows, cols);

	  const MatrixXd &S = Map<MatrixXd>( const_cast<double*>(pS), rows, cols);
	  const MatrixXd I = MatrixXd::Identity(cols,cols);
	  const double n = (S.transpose()*S-I).norm();
	  return n*n*penalty*0.5f;
	}

	void grad(const double *pS, double *pG, const int rows, const int cols)const{

	  assert(pS);
	  assert(pG);
	  assert_ge(rows, cols);

	  const MatrixXd S = Map<MatrixXd>( const_cast<double*>(pS), rows, cols );
	  static MatrixXd g;
	  g = (S*S.transpose()*S)-S;
	  g *= 2*penalty;
	  
	  memcpy(pG, &g(0,0), g.size()*sizeof(double) );
	}

	void hessian(VectorXd &H,const int rw, const int rs)const{
	  // @bug not implemented yet.
	  assert(false);
	}
  };
  typedef boost::shared_ptr<BasisNormalizeEnergy> pBasisNormalizeEnergy;

  class BasisNormalizeColEnergy:public BasisSelRegEnergy{
	
  public:
	double fun(const double *pS, const int rows, const int cols)const{

	  assert(pS);
	  assert_ge(rows, cols);

	  const MatrixXd &S = Map<MatrixXd>( const_cast<double*>(pS), rows, cols);
	  double E = 0.0f;
	  for (int i = 0; i < S.cols(); ++i){
		const double n = S.col(i).dot(S.col(i))-1.0f;
		E += n*n;
	  }
	  return E*penalty*0.5f;
	}

	void grad(const double *pS, double *pG, const int rows, const int cols)const{

	  assert(pS);
	  assert(pG);
	  assert_ge(rows, cols);

	  const MatrixXd S = Map<MatrixXd>( const_cast<double*>(pS), rows, cols );
	  static MatrixXd g;
	  g.resize(S.rows(),S.cols());
	  for (int k = 0; k < S.cols(); ++k){
		g.col(k) = S.col(k)*(S.col(k).dot(S.col(k))-1.0f)*2.0f;
	  }
	  g *= penalty;
	  
	  memcpy(pG, &g(0,0), g.size()*sizeof(double) );
	}

	void hessian(VectorXd &H,const int rw, const int rs)const{
	  // @bug not implemented yet.
	  assert(false);
	}
  };
  
}//end of namespace

#endif /*_BASISSELREGENERGY_H_*/
