#ifndef _WARPINGCUBATURE_H_
#define _WARPINGCUBATURE_H_

#include <assertext.h>
#include <RS2Euler.h>
#include <ComputeBj.h>
#include <MapMA2RS.h>
#include <IO.h>
#include "./cubacode/GreedyCubop.h"
using namespace LSW_WARPING;

namespace LSW_UTILITY{

  typedef std::vector<COMMON::Vec3,Eigen::aligned_allocator<COMMON::Vec3> > VectorV3;
  
  /**
   * @class WarpingCubature base class for warping cubature.
   * 
   */
  class WarpingCubature: public GreedyCubop{

  public:
	WarpingCubature(const MatrixXd &W, /*linearBasis*/
					const MatrixXd &B, /*nonlinearBasis*/
					pVolumetricMesh_const tetmesh){

	  assert(tetmesh != NULL);
	  rs2euler.setTetMesh(tetmesh);

	  const SparseMatrixD &G = rs2euler.get_G();
	  MapMA2RS::computeMapMatPGW(G, W, hatW);

	  const SparseMatrixD &A = rs2euler.get_L();
	  assert_eq(A.cols(),A.rows());
	  assert_eq(A.cols(),B.rows());
	  const MatrixXd T = (B.transpose()*(A*B));
	  P = T.inverse();
	  P = P*(B.transpose()*rs2euler.get_VG_t());
	  this->B = B;
	}
	void setFixedNodes(const vector<int>&fixed_nodes){
	  rs2euler.setFixedNodes(fixed_nodes);
	  rs2euler.precompute();
	}

	const std::vector<int> &getSelPoints()const{
	  return selectedPoints;
	}
	const std::vector<double> &getWeights()const{
	  return weights;
	}
	const VectorXd computeRSDisp(const VectorXd &z){

	  const VectorXd y = hatW*z;
	  VectorXd u;
	  rs2euler.reconstruct(y,u);
	  return u;
	}
	bool saveAsVTK(const string filename)const{
	  
	  VectorV3 points(selectedPoints.size());
	  const pVolumetricMesh_const tetmesh = rs2euler.getTetMesh();
	  for (size_t i = 0; i < selectedPoints.size(); ++i){
	  	const Vec3d v = tetmesh->elementCenter(selectedPoints[i]);
	  	points[i][0] = v[0];
	  	points[i][1] = v[1];
	  	points[i][2] = v[2];
	  }

	  COMMON::VTKWriter<double> writer("point",filename,true);
	  //append point
	  writer.appendPoints(points.begin(),points.end());

	  //append cell
	  COMMON::VTKWriter<double>::IteratorIndex<COMMON::Vec3i> beg(0,0,1);
	  COMMON::VTKWriter<double>::IteratorIndex<COMMON::Vec3i> end(points.size(),0,1);
	  writer.appendCells(beg,end,COMMON::VTKWriter<double>::POINT);

	  std::vector<double> color(points.size());
	  for(size_t i=0;i<points.size();i++)
	  	color[i]=weights[i];
	  writer.appendCustomData("color",color.begin(),color.end());
	  return true;
	}
	void convertTrainingData(const MatrixXd &Z,const MatrixXd &Q,
							 TrainingSet &trZ,VECTOR &trQ){
	  
	  std::vector<VectorXd> Qs;
	  for (int i = 0; i < Z.cols(); ++i){

		const VectorXd &z = Z.col(i);
		const VectorXd &q = Q.col(i);
		if (q.norm() > 1e-5){
		  VECTOR *v = new VECTOR();
		  fromEigenVec(z,*v);
		  trZ.push_back(v);
		  Qs.push_back(q);
		}
	  }

	  assert_gt(Qs.size(),0);
	  VectorXd allQ(Qs[0].size()*Qs.size());
	  for (size_t i = 0; i < Qs.size(); ++i){
		assert_eq(Qs[0].size(),Qs[i].size());
		allQ.segment(Qs[0].size()*i,Qs[0].size()) = Qs[i];
	  }
	  fromEigenVec(allQ,trQ);
	}
	
  protected:
	int numTotalPoints(){
	  pVolumetricMesh_const vol = rs2euler.getTetMesh();
	  if (vol)
		return vol->numElements();
	  return 0;
	}
	int numTotalPoints()const{
	  pVolumetricMesh_const vol = rs2euler.getTetMesh();
	  if (vol)
		return vol->numElements();
	  return 0;
	}
	void handleCubature( std::vector<int>& selPoints,VECTOR& ws,Real relErr){

	  selectedPoints.clear();
	  weights.clear();
	  selectedPoints.reserve(selPoints.size());
	  weights.reserve(selPoints.size());
	  for (size_t i = 0; i < selPoints.size(); ++i){
		if(ws(i) >= 1e-6){
		  selectedPoints.push_back(selPoints[i]);
		  weights.push_back(ws(i));
		}
	  }
	  GreedyCubop::handleCubature(selPoints,ws, relErr);
	}

	void toEigenVec(VECTOR& v, VectorXd &ev)const{
	  ev.resize(v.size());
	  memcpy(&ev[0],v.data(),sizeof(double)*v.size());
	}
	void fromEigenVec(const VectorXd &ev, VECTOR& v)const{
	  v.resizeAndWipe(ev.size());
	  memcpy(v.data(),&ev[0],sizeof(double)*ev.size());
	}

  protected:
	MatrixXd hatW;
	MatrixXd P;
	MatrixXd B;
	RS2Euler rs2euler;
	// results
	std::vector<int> selectedPoints;
	std::vector<double> weights;
  };

}//end of namespace

#endif /*_WARPINGCUBATURE_H_*/
