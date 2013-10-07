#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <UnitTestAssert.h>
#include <DrawCurves.h>
#include <AuxTools.h>
#include <MapMA2RS.h>
#include <MeshVtkIO.h>
#include <MASimulatorAD.h>
#include "ElasticMtlOpt.h"
#include "HarmonicOscillator.h"
#include "DFT.h"
#include <eigen3/Eigen/Dense>
using namespace Eigen;
using namespace LSW_SIM;
using namespace ANI_EDIT;
using namespace UTILITY;
using namespace Eigen;

BOOST_AUTO_TEST_SUITE(ElasticMtlTest)

BOOST_AUTO_TEST_CASE(computeTestUStvk){
  
  const double h = 0.1;
  const int r = 10;
  const int T = 50;
  ElasticMtlOpt mtlOpt(h,r);
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  TEST_ASSERT(mtlOpt.loadVolMesh(data+"/sim-mesh.hinp"));
  MatrixXd U;
  TEST_ASSERT(load(data+"/tempt/UstrCon10.b",U));
  mtlOpt._U.resize(U.rows(),T);
  const int dif = U.cols()/T;
  for (int i = 0; i < T; ++i){
    mtlOpt._U.col(i) = U.col(dif*i);
  }
  cout<< "mtlOptU.norm(): " << mtlOpt._U.norm () << endl; 

  MatrixXd tW,PGW;
  TEST_ASSERT (load(data+"eigen_vectors80.b",tW));
  const MatrixXd W = tW.leftCols(r);
  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(mtlOpt._volMesh,G) );
  MapMA2RS::computeMapMatPGW(G,W,PGW);
  const MatrixXd invPGW = PseudoInverse(PGW);
  cout<< "invPGW.norm() = " << invPGW.norm() << endl;

  MatrixXd Y;
  TEST_ASSERT( load(data+"/tempt/YrsCon10.b",Y) );
  mtlOpt.computeRS();

  cout << "UY.norm: " << mtlOpt._Urs.norm() << endl;
  cout << "Y.norm: " << Y.norm() << endl;
  cout << "Y diff: " << (mtlOpt._Urs-Y).norm() << endl;


  // mtlOpt._zrs.resize(r,T);
  // for (int i = 0; i < T; ++i){
  //   mtlOpt._zrs.col(i) = invPGW*mtlOpt._Urs.col(i);
  // }

  // TEST_ASSERT( load(data+"/tempt/ZrsCon10.b",mtlOpt._zrs) );
  // mtlOpt.computeK();
  // mtlOpt._Wrs = PGW;
  // mtlOpt.decomposeK();
}

BOOST_AUTO_TEST_SUITE_END()
