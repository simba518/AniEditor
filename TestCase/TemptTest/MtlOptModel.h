#ifndef _MTLOPTMODEL_H_
#define _MTLOPTMODEL_H_

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <AuxTools.h>
#include <MatrixIO.h>
#include <RS2Euler.h>
#include <RSCoordComp.h>
#include <DefGradOperator.h>
#include <MapMA2RS.h>
#include <MtlOptEnergyAD.h>
#include <JsonFilePaser.h>
using namespace std;
using namespace UTILITY;
using namespace EIGEN3EXT;
using namespace UTILITY;
using namespace LSW_WARPING;

class MtlOptModel{
  
public:
  MtlOptModel(const string initf);
  void initMtlData(MtlDataModel &model);
  void saveMesh(const MatrixXd &Z,const string fname){
	for (int i = 0; i < Z.cols(); ++i)
	  saveMeshOneZ(Z.col(i),fname+TOSTR(i)+".vtk");
  }
  void saveMeshOneZ(const VectorXd &z,const string fname);
  void saveUc(const string fname)const{
	for (int i = 0; i < uc.size(); ++i)
	  saveUc(fname+TOSTR(i)+".vtk",uc[i],conNodes[i]);
  }
  void saveUc(const string fname,const VectorXd &uc,const vector<int> &nid)const;
  int redDim()const{return lambda.size();}
  void initWarper(JsonFilePaser &jsonf);
  void print()const;

public:
  pTetMesh tetmesh;
  RS2Euler rs2euler;
  pRedRSWarperAD warper;

  MatrixXd W;
  MatrixXd hatW;
  VectorXd lambda;
  vector<int> testModeId;

  int T;
  double h;
  double alphaK;
  double alphaM;

  MatrixXd Z_initial;
  MatrixXd Kz;
  vector<int> Kid;
  vector<int> conFrames;
  vector<vector<int> > conNodes;
  vector<VectorXd> uc;
  double partialConPenalty;
  double fullConPenalty;
};

typedef boost::shared_ptr<MtlOptModel> pMtlOptModel;

#endif /* _MTLOPTMODEL_H_ */
