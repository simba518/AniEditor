#ifndef _MTLOPTENERGY_H_
#define _MTLOPTENERGY_H_

#include <eigen3/Eigen/Sparse>
#include <BaseMtlOptEnergy.h>

namespace LSW_ANI_EDITOR{

  typedef std::vector<Matrix3d,aligned_allocator<Matrix3d> > VMat3d;

  class CtrlForceEnergy:public BaseCtrlForceEnergy,public BaseOptFun{
	
  public:
	void setRedWarper(pRedRSWarperExt warper){
	  _warper = warper;
	}
	int dim()const{
	  return (getT()-_keyframes.size())*reducedDim();
	}
	void init(double *x,const int n){
	  assert_gt(n,0);
	  assert_eq(n,dim());
	  assert_eq(_subZ.size(),n);
	  memcpy(x,&_subZ(0,0),sizeof(double)*n);
	}
	double fun(const double *x);
	void grad(const double *x,double *g);
	void setRlst(const double *x, const double objValue){

	  _objValue = objValue;
	  _subZ.resize(reducedDim(),getT()-_keyframes.size());
	  memcpy(&_subZ(0,0),x,sizeof(double)*_subZ.size());
	  _Z = _subZ;
	  addKeyframes(_Z);
	}
	const MatrixXd &getZ()const{
	  return _Z;
	}

  protected:
	void precompute();
	void addKeyframes(MatrixXd &Z)const;
	void keyframeChanged();
	
  private:
	pRedRSWarperExt _warper;
	Eigen::SparseMatrix<double> _insertKeyZ;
	MatrixXd _subZ;
	MatrixXd _Z;
	VectorXd _Ma;
	VectorXd _Mb;
	VMat3d _H;
  };
  typedef boost::shared_ptr<CtrlForceEnergy> pCtrlForceEnergy;

  class MtlOptEnergy: public BaseMtlOptEnergy{
	
  public:
	void setZ(const MatrixXd &Z);
	double fun(const double *x);
	void grad(const double *x,double *g);

  protected:
	MatrixXd _Z1;
	MatrixXd _Z2v;
	MatrixXd _Z2a;
  };
  typedef boost::shared_ptr<MtlOptEnergy> pMtlOptEnergy;

}//end of namespace

#endif /* _MTLOPTENERGY_H_ */
