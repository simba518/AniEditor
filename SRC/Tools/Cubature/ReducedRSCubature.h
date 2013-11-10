#ifndef _REDUCEDRSCUBATURE_H_
#define _REDUCEDRSCUBATURE_H_

#include "WarpingCubature.h"

namespace LSW_UTILITY{

  /**
   * @class ReducedRSCubature use cubature method to reduce the RS method. Use
   * the code from http://www.cs.cornell.edu/~stevenan/cubature/
   */
  class ReducedRSCubature: public WarpingCubature{
	
  public:
	ReducedRSCubature(const MatrixXd &W,const MatrixXd &B,
					  pVolumetricMesh_const tetmesh):WarpingCubature(W,B,tetmesh){}
	void initTrainingData(const MatrixXd &Z,TrainingSet &trZ,VECTOR &trQ){

	  MatrixXd Q(P.rows(),Z.cols());
	  for (int i = 0; i < Z.cols(); ++i){
		const VectorXd &z = Z.col(i);
		const VectorXd y = hatW*z;
		VectorXd b;
		rs2euler.assemble_b(y,b);
		Q.col(i) = P*b;
	  }
	  convertTrainingData(Z,Q,trZ,trQ);
	}
	const VectorXd computeCubDisp(const VectorXd &z)const{

	  assert_gt(selectedPoints.size(),0);
	  assert_eq(selectedPoints.size(), weights.size());

	  VectorXd q;
	  evalRedDispOfOneTet(selectedPoints[0],z,q);
	  q = q*weights[0];
	  for (size_t i = 1; i < selectedPoints.size(); ++i){
		VectorXd tq;
		evalRedDispOfOneTet(selectedPoints[i],z,tq);
		q = q + weights[i]*tq;
	  }
	  return q;
	}
	
  protected:
	void evalPointForceDensity( int pointId, VECTOR& _z, VECTOR& _gOut ){

	  VectorXd z,q;
	  toEigenVec(_z,z);
	  evalRedDispOfOneTet(pointId,z,q);
	  fromEigenVec(q, _gOut);
	}
	void evalRedDispOfOneTet(int tetId, const VectorXd &z, VectorXd &q)const{

	  assert_in(tetId,0,numTotalPoints()-1);
	  assert_eq(numTotalPoints()*9, P.cols());
	  assert_eq(hatW.cols(),z.size());
	  assert_eq(hatW.rows(),numTotalPoints()*9);
	  
	  const VectorXd yj = hatW.block(tetId*9,0,9,hatW.cols())*z;
	  VectorXd bj;
	  const vector<double> &Sqrt_V = rs2euler.get_Sqrt_V();
	  ComputeBj::compute(yj,Sqrt_V[tetId],bj);
	  q = P.block(0,tetId*9,P.rows(),9)*bj;
	}
  };
  
}//end of namespace

#endif /*_REDUCEDRSCUBATURE_H_*/
