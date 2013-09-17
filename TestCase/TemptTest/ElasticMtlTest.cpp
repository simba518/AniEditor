#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <DrawCurves.h>
#include <MapMA2RS.h>
#include <MeshVtkIO.h>
#include <VolObjMesh.h>
#include <MASimulatorAD.h>
#include "ElasticMtlOpt.h"
#include "HarmonicOscillator.h"
using namespace LSW_SIM;
using namespace ANI_EDIT;
using namespace UTILITY;


BOOST_AUTO_TEST_SUITE(ElasticMtlTest)

BOOST_AUTO_TEST_CASE(symIndexTest){

  MatrixXd M(3,3);
  M << 1,2,3,  2,4,5, 3,5,6;
  ASSERT_EQ(M,M.transpose());
  
  VectorXd v(6);
  v << 1,2,4,3,5,6;
  for (int i = 0; i < 3; ++i){
    for (int j = 0; j < 3; ++j){
	  const int n = ElasticMtlOpt::symIndex(i,j);
	  ASSERT_GE(n,0);
	  ASSERT_LT(n,v.size());
	  ASSERT_EQ(v[n],M(i,j));
	}
  }  
}

BOOST_AUTO_TEST_CASE(produceSymetricMatTest){

  vector<CasADi::SX> s;
  CasADi::SXMatrix K;
  const int dim = 2;
  ElasticMtlOpt::produceSymetricMat("x",dim,s,K);
  ASSERT_EQ(K.elem(0,0).toString(),std::string("x_0"));
  ASSERT_EQ(K.elem(0,1).toString(),std::string("x_1"));
  ASSERT_EQ(K.elem(1,0).toString(),std::string("x_1"));
  ASSERT_EQ(K.elem(1,1).toString(),std::string("x_2"));

  const int dim1 = 1;
  ElasticMtlOpt::produceSymetricMat("x",dim1,s,K);
  ASSERT_EQ(K.size1(),1);
  ASSERT_EQ(K.size2(),1);
  ASSERT_EQ(K.elem(0,0).toString(),std::string("x_0"));
}

BOOST_AUTO_TEST_CASE(checkPCAuseOneMode){

  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  const string hatWstr = data+"/tempt/PGW.b";
  MatrixXd hatW;
  TEST_ASSERT( load(hatWstr,hatW) );

  const int T = 100;
  const int mid = 6;
  
  ElasticMtlOpt mtlotp(0.1,1);
  mtlotp._Urs.resize(hatW.rows(),T);
  VectorXd z(T);
  for (int i = 0; i < T; ++i){
	z[i] = 0.02f*i;
	mtlotp._Urs.col(i) = hatW.col(mid)*z[i];
  }
  mtlotp.PCAonRS();
  cout << (mtlotp._zrs.row(0).transpose()-z*hatW.col(mid).norm()).norm() << endl;
  cout << (hatW.col(mid)/hatW.col(mid).norm()-mtlotp._Wrs.col(0)).norm() << endl;
}

BOOST_AUTO_TEST_CASE(checkPCAuseTwoModes){

  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  const string hatWstr = data+"/tempt/PGW.b";
  MatrixXd hatW;
  TEST_ASSERT( load(hatWstr,hatW) );

  const int T = 100;
  const int mid0 = 6;
  const int mid1 = 7;
  
  ElasticMtlOpt mtlotp(0.1,2);
  mtlotp._Urs.resize(hatW.rows(),T);
  MatrixXd z(2,T);
  for (int i = 0; i < T; ++i){
	z(0,i) = 0.02f*i;
	z(1,i) = 0.005f*i;
	mtlotp._Urs.col(i) = hatW.col(mid0)*z(0,i)+hatW.col(mid1)*z(1,i);
  }
  mtlotp.PCAonRS();

  cout<<mtlotp._Wrs.col(0).transpose()*mtlotp._Wrs.col(1)<<endl;
  cout<<hatW.col(mid0).normalized().transpose()*hatW.col(mid1).normalized()<<endl;

  cout<<(mtlotp._zrs.row(0)-z.row(0)*hatW.col(mid0).norm()).norm() << endl;
  cout<<(mtlotp._zrs.row(1)-z.row(1)*hatW.col(mid1).norm()).norm() << endl;

  cout<<(hatW.col(mid0).normalized()-mtlotp._Wrs.col(0)).norm()<<endl;
  cout<<(hatW.col(mid1).normalized()-mtlotp._Wrs.col(1)).norm()<<endl;
}

BOOST_AUTO_TEST_CASE(computeKuseRSSimulation){

  LSW_ANI_EDITOR::MASimulatorAD simulator;
  const double h = 0.1f;
  const int r = 5;
  VectorXd lambda(r);
  for (int i = 0; i < r; ++i){
	lambda[i] = 0.5f*(i+1.0f);
  }
  simulator.setTimeStep(h);
  simulator.setEigenValues(lambda);
  simulator.setStiffnessDamping(0.0f);
  simulator.setMassDamping(0.0f);
  const VectorXd v0 = VectorXd::Zero(r);
  VectorXd z0(r);
  for (int i = 0; i < r; ++i){
	z0[i] = 1.0f/(i+1.0f);
  }
  simulator.setIntialStatus(v0,z0);

  const int T = 50;
  ElasticMtlOpt mtlopt(h,r);
  mtlopt._zrs.resize(r,T);
  for (int f = 0; f < T; ++f){
	mtlopt._zrs.col(f) = simulator.getEigenZ();
    simulator.forward(VectorXd::Zero(r));
  }
  mtlopt.computeK();

  // const VectorXd z = mtlopt._zrs.row(0).transpose();
  // TEST_ASSERT( PythonScriptDraw2DCurves::write("t.py",z,h) );  

  // HarmonicOscillator<double> fz(mtlopt._K(0,0),mtlopt._D[0],z0[0],v0[0]);
  // const VectorXd z2 = fz.generateSequence<VectorXd>(0.0f,h,z.size());
  // TEST_ASSERT( PythonScriptDraw2DCurves::write("t2.py",z2,h) );
}

BOOST_AUTO_TEST_CASE(computeTestUrs){
  
  ElasticMtlOpt mtlOpt(0.1f,30);
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  TEST_ASSERT(mtlOpt.loadVolMesh(data+"/sim-mesh.hinp"));
  TEST_ASSERT(mtlOpt.loadAniSeq(data+"/tempt/UstrCon10.b"));
  cout<< "mtlOptU.norm(): " << mtlOpt._U.norm () << endl;
  mtlOpt.compute();
}

BOOST_AUTO_TEST_CASE(computeTestUStvk){
  
  const double h = 0.1;
  const int r = 10;
  const int T = 100;
  ElasticMtlOpt mtlOpt(h,r);
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  TEST_ASSERT(mtlOpt.loadVolMesh(data+"/sim-mesh.hinp"));
  MatrixXd U;
  TEST_ASSERT(load(data+"/tempt/Ustvk_nocon_len1602.b",U));
  mtlOpt._U.resize(U.rows(),T);
  
  const int dif = U.cols()/T;
  for (int i = 0; i < T; ++i){
    mtlOpt._U.col(i) = U.col(dif*i);
  }
  
  cout<< "mtlOptU.norm(): " << mtlOpt._U.norm () << endl; 

  MatrixXd tW,PGW;
  TEST_ASSERT (load(data+"eigen_vectors80.b",tW));
  const MatrixXd W = tW.leftCols(r);
  mtlOpt.computeRS();

  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(mtlOpt._volMesh,G) );
  MapMA2RS::computeMapMatPGW(G,W,PGW);
  const MatrixXd invPGW = PseudoInverse(PGW);
  cout << "(invPGW*PGW-I).norm() = " << (invPGW*PGW-MatrixXd::Identity(r,r)).norm() << endl;

  cout<< "invPGW.norm() = " << invPGW.norm() << endl;

  mtlOpt._zrs.resize(r,T);
  for (int i = 0; i < T; ++i){
    mtlOpt._zrs.col(i) = invPGW*mtlOpt._Urs.col(i);
  }

  mtlOpt.computeK();
  mtlOpt._Wrs = PGW;
  mtlOpt.decomposeK();

  // save animation.
  VolObjMesh volobj;
  volobj.loadObjMesh(data+"beam.obj");
  volobj.loadVolMesh(data+"sim-mesh.hinp");
  volobj.loadInterpWeights(data+"interp-weights.txt");

  // vector<int> fixed_nodes;
  // TEST_ASSERT( loadVec(data+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
  
  RS2Euler rs2euler;
  rs2euler.setTetMesh(volobj.getVolMesh());
  rs2euler.fixBaryCenter();
  // rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();

  for (int i = 0; i < T; ++i){
	const VectorXd u = mtlOpt._U.col(i);
  	TEST_ASSERT ( volobj.interpolate(&u[0]) );
  	const string filename = data+"tempt/stvk_obj_"+MeshVtkIO::int2str(i)+".vtk";
  	TEST_ASSERT (MeshVtkIO::write(volobj.getObjMesh(),filename));
  }

}

BOOST_AUTO_TEST_CASE(produceUrs){

  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  const string eval = data+"/eigen_values80.b";
  const string evect = data+"/eigen_vectors80.b";
  const int r = 10;
  VectorXd tlambda;
  TEST_ASSERT( load(eval,tlambda) );
  const VectorXd lambda = tlambda.segment(0,r);
  MatrixXd tW;
  TEST_ASSERT( load(evect,tW) );
  const MatrixXd W = tW.leftCols(r);

  const string z0str = data+"/tempt/z135.b";
  VectorXd tz0;
  TEST_ASSERT( load(z0str,tz0) );
  const VectorXd z0 = tz0.segment(0,r);
  
  LSW_ANI_EDITOR::MASimulatorAD simulator;
  const double h = 0.1f;
  simulator.setTimeStep(h);
  simulator.setEigenValues(lambda);
  simulator.setStiffnessDamping(0.0f);
  simulator.setMassDamping(0.0f);
  const VectorXd v0 = VectorXd::Zero(r);
  simulator.setIntialStatus(v0,z0);

  const int T = 50;
  ElasticMtlOpt mtlopt(h,r);
  mtlopt._zrs.resize(r,T);
  for (int f = 0; f < T; ++f){
  	mtlopt._zrs.col(f) = simulator.getEigenZ();
    simulator.forward(VectorXd::Zero(r));
  }

  cout<< endl << lambda << endl;

  // mtlopt.computeK();
  
  // save animation.
  VolObjMesh volobj;
  volobj.loadObjMesh(data+"beam.obj");
  volobj.loadVolMesh(data+"sim-mesh.hinp");
  volobj.loadInterpWeights(data+"interp-weights.txt");

  vector<int> fixed_nodes;
  TEST_ASSERT( loadVec(data+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
  
  RS2Euler rs2euler;
  rs2euler.setTetMesh(volobj.getVolMesh());
  // rs2euler.fixBaryCenter();
  rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();

  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(volobj.getVolMesh(),G) );
  VectorXd y,u;
  vector<VectorXd> U;
  for (int i = 0; i < T; ++i){
  	const VectorXd p = W*mtlopt._zrs.col(i);
  	RSCoordComp::constructWithWarp(G,p,y);
  	TEST_ASSERT ( rs2euler.reconstruct(y,u) );
  	TEST_ASSERT ( volobj.interpolate(&u[0]) );
  	const string filename = data+"tempt/obj_"+MeshVtkIO::int2str(i)+".vtk";
  	TEST_ASSERT (MeshVtkIO::write(volobj.getObjMesh(),filename));
	U.push_back(u);
  }
  
  const string Ustr = data+"/tempt/UstrCon"+MeshVtkIO::int2str(r)+".b";
  TEST_ASSERT( write(Ustr,U));
  
  TEST_ASSERT( write(data+"/tempt/u0.b",U[0] ));
}

BOOST_AUTO_TEST_SUITE_END()

// BOOST_AUTO_TEST_CASE(produceCorrectData){
  
//   const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
//   const string volstr = data+"/sim-mesh.hinp";
//   pVolumetricMesh volMesh;
//   volMesh=pVolumetricMesh(VolumetricMeshLoader::load(volstr.c_str()));
//   TEST_ASSERT(volMesh != NULL);

//   const string wstr = data+"eigen_vectors80.b";
//   MatrixXd W_t;
//   TEST_ASSERT( load(wstr,W_t));
//   const int r = 20;
//   const MatrixXd W = W_t.leftCols(r);

//   SparseMatrix<double> G;
//   TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(volMesh,G) );

//   MatrixXd PGW;
//   LSW_WARPING::MapMA2RS::computeMapMatPGW(G,W,PGW);
//   const string pgw = data+"/tempt/PGW.b";
//   TEST_ASSERT( write(pgw,PGW) );
  
//   const MatrixXd invPGW = PseudoInverse(PGW);
//   const string invpgw = data+"/tempt/invPGW.b";
//   TEST_ASSERT( write(invpgw,invPGW));
// }

// BOOST_AUTO_TEST_CASE(computeKuseRSSimulation){

//   LSW_ANI_EDITOR::MASimulatorAD simulator;
//   const double h = 0.1f;
//   const int r = 10;
//   VectorXd lambda(r);
//   for (int i = 0; i < r; ++i){
// 	lambda[i] = 0.5f*(i+1.0f);
//   }
//   simulator.setTimeStep(h);
//   simulator.setEigenValues(lambda);
//   simulator.setStiffnessDamping(0.0f);
//   simulator.setMassDamping(0.0f);
//   const VectorXd v0 = VectorXd::Zero(r);
//   VectorXd z0(r);
//   for (int i = 0; i < r; ++i){
// 	z0[i] = 10.0f/(i+1.0f);
//   }
//   simulator.setIntialStatus(v0,z0);

//   const int T = 100;
//   ElasticMtlOpt mtlopt(h*10.0f,r);
//   mtlopt._zrs.resize(r,T);
//   for (int f = 0; f < T; ++f){
// 	mtlopt._zrs.col(f) = simulator.getEigenZ();
//     simulator.forward(VectorXd::Zero(r));
//   }
//   mtlopt.computeK();
// }

// BOOST_AUTO_TEST_CASE(computeTest){
  
//   ElasticMtlOpt mtlOpt(0.1f,20);
//   const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
//   TEST_ASSERT(mtlOpt.loadVolMesh(data+"/sim-mesh.hinp"));
//   TEST_ASSERT(mtlOpt.loadAniSeq(data+"mtlOptU.b"));
//   const MatrixXd U = mtlOpt._U.leftCols(50);
//   mtlOpt._U = U;
//   cout<< "mtlOptU.norm(): " << mtlOpt._U.norm () << endl;
//   mtlOpt.compute();
// }
