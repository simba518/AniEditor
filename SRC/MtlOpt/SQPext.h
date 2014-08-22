#ifndef HJ_SQP_H_
#define HJ_SQP_H_

#include <linear_solver.h>

namespace MTLOPT {

  class SQPfun{
	
  public:
	virtual size_t dim(void) const = 0;

	virtual int prepare_for_new_x(const double *x) = 0;

	virtual int val(const double *x, double &v) const = 0;

	virtual int gra(const double *x, double *g) = 0;

	//NOTE: add to hes
	//! @param h,ptr,idx == 0 means query nnz
	//! @param h == 0, ptr,idx != 0 means query patten
	//! @param h,ptr,idx != 0 means accumulate
	//! @param format for query pattern: 1 csc, 2 pair of row, col
	virtual int hes(const double *x, size_t &nnz, size_t &format, double *h,
					int *ptr, int *idx, double alpha = 1) = 0;

  };

  class SQPext{

  public:
    SQPext():f_(0){}

	void clear_hes(){
	  H_.resize(0, 0, 0);
	}

    int set_f(SQPfun *f);

    int solve(double *x, boost::property_tree::ptree &pt);

  protected:
    SQPfun *f_;

    std::unique_ptr<linear_solver> slv_;

    hj::sparse::csc<double, int32_t> H_;

    zjucad::matrix::matrix<double> g_, s_, tmp_;
  };

}

#endif
