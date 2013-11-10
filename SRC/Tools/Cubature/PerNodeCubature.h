#ifndef _PERNODECUBATURE_H_
#define _PERNODECUBATURE_H_

#include "WarpingCubature.h"

namespace LSW_UTILITY{

  /**
   * @class PerNodeCubature use cubature to compute the sample points and
   * weights for the warping of each node.
   * @see FittingRSCurve.pdf
   */
  class PerNodeCubature: public WarpingCubature{

  public:	
	PerNodeCubature(const MatrixXd &W, const MatrixXd &B, 
					pVolumetricMesh_const tetmesh):
	  WarpingCubature(W,B,tetmesh){
	  node_to_cubature = 0;
	}
	void setNodeToCubature(const int node){
	  assert_in(node*3,0,hatW.rows()-3);
	  this->node_to_cubature = node;
	}
	void initTrainingData(const MatrixXd &Z,TrainingSet &trZ,VECTOR &trU3){

	  MatrixXd U(3,Z.cols());
	  for (int i = 0; i < Z.cols(); ++i){
		const VectorXd &z = Z.col(i);
		const VectorXd y = hatW*z;
		VectorXd b;
		rs2euler.assemble_b(y,b);
		U.col(i) = B.block(nodeToCubature()*3,0,3,B.cols())*(P*b);
	  }
	  convertTrainingData(Z,U,trZ,trU3);
	}
	int nodeToCubature()const{
	  return node_to_cubature;
	}
	
  protected:
	void evalPointForceDensity( int pointId, VECTOR& _z, VECTOR& _gOut ){

	  VectorXd z;
	  Vector3d u3;
	  toEigenVec(_z,z);
	  evalDispOfOneTet(pointId,z,u3);
	  assert_eq(_gOut.size(),3);
	  _gOut(0) = u3[0];
	  _gOut(1) = u3[1];
	  _gOut(2) = u3[2];
	}
	void evalDispOfOneTet(int tetId, const VectorXd &z, Vector3d &u3)const{

	  assert_in(tetId,0,numTotalPoints()-1);
	  assert_eq(numTotalPoints()*9, P.cols());
	  assert_eq(hatW.cols(),z.size());
	  assert_eq(hatW.rows(),numTotalPoints()*9);
	  
	  const VectorXd yj = hatW.block(tetId*9,0,9,hatW.cols())*z;
	  VectorXd bj;
	  const vector<double> &Sqrt_V = rs2euler.get_Sqrt_V();
	  ComputeBj::compute(yj,Sqrt_V[tetId],bj);
	  const VectorXd q = P.block(0,tetId*9,P.rows(),9)*bj;
	  u3 = B.block(nodeToCubature()*3,0,3,B.cols())*q;
	}

  private:
	int node_to_cubature;
  };
  
}//end of namespace

#endif /*_PERNODECUBATURE_H_*/
