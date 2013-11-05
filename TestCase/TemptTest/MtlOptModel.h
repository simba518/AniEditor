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
#include <MeshVtkIO.h>
#include <limits>
#include <MtlOptEnergyAD.h>
using namespace std;
using namespace EIGEN3EXT;
using namespace UTILITY;
using namespace LSW_WARPING;

typedef struct _MtlOptModel{

  _MtlOptModel(const string initf);
  bool loadLambda(const string initf);
  void initVolObj(const string initf);
  void produceSimRlst(const bool genKeyZ=true);
  void extrangeKeyframes();
  void initMtlOpt(RedSpaceTimeEnergyAD &ad)const;
  void initMtlData(MtlDataModel &model);
  void initSolver(const SXMatrix &E, const VSX &x);
  void solve();
  VectorXd getOutput();
  void getZfromSolver(MatrixXd &Z);
  void saveRlst(const string dir);
  void saveMesh(const MatrixXd &Z,const string dir,const string fname);
  void saveMeshOneZ(const VectorXd &z,const string fname);
  void computeEnergy(const VectorXd &X);
  int redDim()const;
  static MatrixXd assembleFullZ(const VectorXd&subZ,const VectorXd&keyZ,const VectorXi&Kid, const int r);

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

  CasADi::SXFunction fun;
  CasADi::IpoptSolver solver;
  VSX allVars;

}MtlOptModel;

#endif /* _MTLOPTMODEL_H_ */
