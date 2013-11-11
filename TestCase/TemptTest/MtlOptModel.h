#ifndef _MTLOPTMODEL_H_
#define _MTLOPTMODEL_H_

#include <iostream>
#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <interfaces/ipopt/ipopt_solver.hpp>
#include <AuxTools.h>
#include <MatrixIO.h>
#include <MASimulatorAD.h>
#include <TetMeshEmbeding.h>
#include <RS2Euler.h>
#include <RSCoordComp.h>
#include <DefGradOperator.h>
#include <MapMA2RS.h>
#include <limits>
#include <MtlOptEnergyAD.h>
#include <JsonFilePaser.h>
using namespace UTILITY;
using namespace std;
using namespace EIGEN3EXT;
using namespace UTILITY;
using namespace LSW_WARPING;

typedef struct _MtlOptModel{

  _MtlOptModel(const string initf);
  bool loadLambda(const string initf);
  void initVolObj(const string initf);
  void produceSimRlst(const bool genKeyZ=true);
  void initMtlData(MtlDataModel &model);
  void saveMesh(const MatrixXd &Z,const string fname);
  void saveMeshOneZ(const VectorXd &z,const string fname);
  void saveUc(const string fname)const;
  void saveUc(const string fname,const VectorXd &uc,const vector<int> &nid)const;
  int redDim()const;
  static MatrixXd assembleFullZ(const VectorXd&subZ,const VectorXd&keyZ,const VectorXi&Kid, const int r);
  void initWarper(JsonFilePaser &jsonf);
  void print()const;

  MatrixXd W;
  SparseMatrix<double> G;
  MatrixXd PGW;
  TetMeshEmbeding volobj;
  RS2Euler rs2euler;

  int T;
  double h;
  double alphaK;
  double alphaM;
  VectorXd lambda;
  VectorXi testModeId;
  VectorXd z0;
  MatrixXd Z;
  MatrixXd CorrectZ;

  MatrixXd Kz;
  VectorXi Kid;
  vector<int> conFrames;
  vector<vector<int> > conNodes;
  vector<VectorXd> uc;
  double penaltyCon;
  pRedRSWarperAD warper;

}MtlOptModel;
typedef boost::shared_ptr<MtlOptModel> pMtlOptModel;

#endif /* _MTLOPTMODEL_H_ */
