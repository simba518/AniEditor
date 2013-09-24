#ifndef _ELASTICMTLOPT_H_
#define _ELASTICMTLOPT_H_

#include <string>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include <CASADITools.h>
#include <MatrixIO.h>
#include <MatrixTools.h>
#include <DefGradOperator.h>
#include <RSCoordComp.h>
#include <RS2Euler.h>
#include <volumetricMesh.h>
#include <volumetricMeshLoader.h>
using namespace std;
using namespace Eigen;
using namespace EIGEN3EXT;
using namespace CASADI;  
using namespace LSW_WARPING;
using namespace UTILITY;  

namespace ANI_EDIT{

  class ElasticMtlOpt{
  
  public:
	ElasticMtlOpt(double h=0.1f,int PCADim=80, double ak=0.0f,double am=0.0f):
	  _h(h),_alphaK(ak),_alphaM(am),_PCADim(PCADim){}
	bool loadVolMesh(const string filename){
	  _volMesh = pVolumetricMesh(VolumetricMeshLoader::load(filename.c_str()));
	  return (_volMesh != NULL);
	}
	bool loadAniSeq(const string filename){
	  return EIGEN3EXT::load(filename,_U);
	}
	void compute(){
	  computeRS();
	  PCAonRS();
	  computeK();
	  decomposeK();
	}
  
	void computeRS();
	void PCAonRS();
	void computeK();
	void decomposeK();
	static void produceSymetricMat(const string name,const int dim,vector<CasADi::SX> &s,CasADi::SXMatrix &SM);
	static int symIndex(int r,int c);

	int getT()const{
	  return _U.cols();
	}
  
  public:
	double _h;
	double _alphaK;
	double _alphaM;
	int _PCADim;
	pVolumetricMesh _volMesh;
	MatrixXd _U;

	MatrixXd _Urs;
	MatrixXd _Wrs;
	MatrixXd _zrs;
	MatrixXd _K;
	MatrixXd _D;
	MatrixXd _B;
	VectorXd _Lambda;
	MatrixXd _W;

	CasADi::SXFunction _Efun;
  };
}

#endif /* _ELASTICMTLOPT_H_ */
