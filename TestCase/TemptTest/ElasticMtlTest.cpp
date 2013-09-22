#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <UnitTestAssert.h>
#include <DrawCurves.h>
#include <MapMA2RS.h>
#include <MeshVtkIO.h>
#include <VolObjMesh.h>
#include <MASimulatorAD.h>
#include "ElasticMtlOpt.h"
#include "HarmonicOscillator.h"
#include "DFT.h"
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
  const double h = 0.5f;
  const int r = 1;
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

  PythonScriptDraw2DCurves<VectorXd> curves;
  const VectorXd z = mtlopt._zrs.row(0).transpose();
  curves.add("z_old",z,h,0.0f,"--");

  HarmonicOscillator<double> fz(mtlopt._K(0,0),mtlopt._D[0],z0[0],v0[0]);
  const VectorXd z2 = fz.generateSequence<VectorXd>(0.0f,h,z.size());
  curves.add("z_new",z2,h,0.0f,"*");
  TEST_ASSERT( curves.write("a.py") );
  
  cout << mtlopt._K(0,0) << endl << endl;
  VectorXd outreal,outimag;
  computeDFT(z,outreal,outimag);
  curves.clear();
  curves.add("real",outreal,h,0,"o");
  curves.add("imag",outimag,h,0,"*");
  TEST_ASSERT( curves.write("b.py") );

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

BOOST_AUTO_TEST_CASE(produceUrs){

  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  // const string eval = data+"/eigen_values80.b";
  // const string evect = data+"/eigen_vectors80.b";
  const string eval = data+"/eigen_values_nocon80.b";
  const string evect = data+"/eigen_vectors_nocon80.b";
  const int r = 10;
  VectorXd tlambda;
  TEST_ASSERT( load(eval,tlambda) );
  const VectorXd lambda = tlambda.segment(0,r);
  MatrixXd tW;
  TEST_ASSERT( load(evect,tW) );
  const MatrixXd W = tW.leftCols(r);

  // const string z0str = data+"/tempt/z135.b";
  const string z0str = data+"/tempt/z84-nocon.b";
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

  const int T = 100;
  ElasticMtlOpt mtlopt(h,r);
  mtlopt._zrs.resize(r,T);
  for (int f = 0; f < T; ++f){
  	mtlopt._zrs.col(f) = simulator.getEigenZ();
    simulator.forward(VectorXd::Zero(r));
  }
  
  // save animation.
  VolObjMesh volobj;
  volobj.loadObjMesh(data+"beam.obj");
  volobj.loadVolMesh(data+"sim-mesh.hinp");
  volobj.loadInterpWeights(data+"interp-weights.txt");

  vector<int> fixed_nodes;
  TEST_ASSERT( loadVec(data+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
  
  RS2Euler rs2euler;
  rs2euler.setTetMesh(volobj.getVolMesh());
  // rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.fixBaryCenter();
  rs2euler.precompute();

  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(volobj.getVolMesh(),G) );
  MatrixXd PGW;
  MapMA2RS::computeMapMatPGW(G,W,PGW);
  const MatrixXd invPGW = PseudoInverse(PGW);
  cout << "(invPGW*PGW-I).norm() = " << (invPGW*PGW-MatrixXd::Identity(r,r)).norm() << endl;
  cout<< "invPGW.norm() = " << invPGW.norm() << endl;

  VectorXd y,u;
  vector<VectorXd> U,Z,Y;
  for (int i = 0; i < T; ++i){

  	const VectorXd p = W*mtlopt._zrs.col(i);
  	RSCoordComp::constructWithWarp(G,p,y);
	Y.push_back(y);

	const VectorXd z = invPGW*y;
	cout<< "d: " << (z-mtlopt._zrs.col(i)).norm()/mtlopt._zrs.col(i).norm() << endl;
	Z.push_back(z);

  	TEST_ASSERT ( rs2euler.reconstruct(y,u) );
  	TEST_ASSERT ( volobj.interpolate(&u[0]) );
  	const string filename = data+"tempt/obj_"+MeshVtkIO::int2str(i)+".vtk";
  	TEST_ASSERT (MeshVtkIO::write(volobj.getObjMesh(),filename));
	U.push_back(u);
  }
  
  const string Ustr = data+"/tempt/Urs"+MeshVtkIO::int2str(r)+".b";
  TEST_ASSERT( write(Ustr,U));

  const string Zstr = data+"/tempt/Zrs"+MeshVtkIO::int2str(r)+".b";
  TEST_ASSERT( write(Zstr,Z));

  const string Ystr = data+"/tempt/Yrs"+MeshVtkIO::int2str(r)+".b";
  TEST_ASSERT( write(Ystr,Y));
  
  TEST_ASSERT( write(data+"/tempt/u0.b",U[0] ));
}

BOOST_AUTO_TEST_CASE(checkRS){
  
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  VolObjMesh volobj;
  volobj.loadVolMesh(data+"sim-mesh.hinp");
  vector<int> fixed_nodes;
  TEST_ASSERT( loadVec(data+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
  
  RS2Euler rs2euler;
  rs2euler.setTetMesh(volobj.getVolMesh());
  rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();

  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(volobj.getVolMesh(),G) );

  const int r = 10;
  MatrixXd tW,PGW;
  TEST_ASSERT ( load(data+"eigen_vectors80.b",tW) );
  const MatrixXd W = tW.leftCols(r);
  MapMA2RS::computeMapMatPGW(G,W,PGW);
  const MatrixXd invPGW = PseudoInverse(PGW);

  VectorXd u1,u2;
  vector<VectorXd> Y1,Y2;
  TEST_ASSERT( load(data+"/tempt/YrsCon10.b",Y1) );

  const int T = Y1.size();
  MatrixXd Z1(r,T),Z2(r,T);
  for (size_t i = 0; i < Y1.size(); ++i){

	const VectorXd &y1 = Y1[i];
	TEST_ASSERT ( rs2euler.reconstruct(y1,u1) );
	VectorXd y2;
	LSW_WARPING::RSCoordComp::constructWithoutWarp(G,u1,y2);
	TEST_ASSERT ( rs2euler.reconstruct(y2,u2) );
	Y2.push_back(y2);
	Z1.col(i) = invPGW*y1;
	Z2.col(i) = invPGW*y2;

	cout<< "i = " << i;
	cout<< "\t yd :"<<(Y2[i]-Y1[i]).norm()/Y1[i].norm();
	cout<< "\t zd :"<<(Z2.col(i)-Z1.col(i)).norm()/Z1.col(i).norm()<<endl;
  }

  const int start = 20;
  const double h = 0.1f;
  ElasticMtlOpt mtlOpt(h,r);
  mtlOpt._zrs.resize(r,T-start);
  for (int i = 0; i < mtlOpt._zrs.cols(); ++i){

    mtlOpt._zrs.col(i) = Z2.col(i+start);
	cout<< "i = " << i;
	cout<< "\t yd :"<<(Y2[i+start]-Y1[i+start]).norm()/Y1[i+start].norm();
	cout<< "\t zd :"<<(Z2.col(i+start)-Z1.col(i+start)).norm()/Z1.col(i+start).norm()<<endl;
  }

  // mtlOpt._zrs = Z2.block(2,20,8,70);
  mtlOpt.computeK();

  VectorXd eigenval;
  TEST_ASSERT ( load(data+"eigen_values80.b",eigenval) );
  cout<< "correct eigen values:\n" << eigenval.segment(0,r) << endl;;

  mtlOpt._Wrs = PGW;
  mtlOpt.decomposeK();
 

  PythonScriptDraw2DCurves<VectorXd> curves;
  for (int i = 0; i < 3; ++i){
	const string index = boost::lexical_cast<string>(i);
	TEST_ASSERT( curves.add(string("old")+index,(VectorXd)(Z1.row(i)),h,0,"o") );
	TEST_ASSERT( curves.add(string("new")+index,(VectorXd)(Z2.row(i)),h,0) );
  }
  TEST_ASSERT ( curves.write("t1.py") );
  
  curves.clear();
  for (int i = 3; i < 6; ++i){
	const string index = boost::lexical_cast<string>(i);
	TEST_ASSERT( curves.add(string("old")+index,(VectorXd)(Z1.row(i)),h,0,"o") );
	TEST_ASSERT( curves.add(string("new")+index,(VectorXd)(Z2.row(i)),h,0) );
  }
  TEST_ASSERT ( curves.write("t2.py") );

  curves.clear();
  for (int i = 6; i < 10; ++i){
	const string index = boost::lexical_cast<string>(i);
	TEST_ASSERT( curves.add(string("old")+index,(VectorXd)(Z1.row(i)),h,0,"o") );
	TEST_ASSERT( curves.add(string("new")+index,(VectorXd)(Z2.row(i)),h,0) );
  }
  TEST_ASSERT ( curves.write("t3.py") );
}

BOOST_AUTO_TEST_CASE(checkRS2){
  
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  VolObjMesh volobj;
  volobj.loadVolMesh(data+"sim-mesh.hinp");

  RS2Euler rs2euler;
  rs2euler.setTetMesh(volobj.getVolMesh());
  rs2euler.fixBaryCenter();
  rs2euler.precompute();

  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(volobj.getVolMesh(),G) );

  const int r = 10;
  MatrixXd tW,PGW;
  TEST_ASSERT ( load(data+"eigen_vectors_nocon80.b",tW) );
  const MatrixXd W = tW.leftCols(r);
  MapMA2RS::computeMapMatPGW(G,W,PGW);
  const MatrixXd invPGW = PseudoInverse(PGW);

  VectorXd u1,u2;
  vector<VectorXd> Y1,Y2;
  TEST_ASSERT( load(data+"/tempt/Yrs10.b",Y1) );

  const int T = Y1.size();
  MatrixXd Z1(r,T),Z2(r,T);
  for (size_t i = 0; i < Y1.size(); ++i){

	const VectorXd &y1 = Y1[i];
	TEST_ASSERT ( rs2euler.reconstruct(y1,u1) );
	VectorXd y2;
	LSW_WARPING::RSCoordComp::constructWithoutWarp(G,u1,y2);
	TEST_ASSERT ( rs2euler.reconstruct(y2,u2) );
	Y2.push_back(y2);
	Z1.col(i) = invPGW*y1;
	Z2.col(i) = invPGW*y2;
  }

  const int start = 20;
  const double h = 0.1f;
  ElasticMtlOpt mtlOpt(h,r);
  mtlOpt._zrs.resize(r,T-start);
  for (int i = 0; i < mtlOpt._zrs.cols(); ++i){

    mtlOpt._zrs.col(i) = Z2.col(i+start);
	cout<< "i = " << i;
	cout<< "\t yd :"<<(Y2[i+start]-Y1[i+start]).norm()/Y1[i+start].norm();
	cout<< "\t zd :"<<(Z2.col(i+start)-Z1.col(i+start)).norm()/Z1.col(i+start).norm()<<endl;
  }
  mtlOpt.computeK();
  mtlOpt._Wrs = PGW;
  mtlOpt.decomposeK();
  
  VectorXd eigenval;
  TEST_ASSERT ( load(data+"eigen_values_nocon80.b",eigenval) );
  cout<< "correct eigen values:\n" << eigenval.segment(0,r) << endl;;

  PythonScriptDraw2DCurves<VectorXd> curves;
  for (int i = 0; i < 3; ++i){
	const string index = boost::lexical_cast<string>(i);
	TEST_ASSERT( curves.add(string("old")+index,(VectorXd)(Z1.row(i)),h,0,"o") );
	TEST_ASSERT( curves.add(string("new")+index,(VectorXd)(Z2.row(i)),h,0) );
  }
  TEST_ASSERT ( curves.write("tt1.py") );
  
  curves.clear();
  for (int i = 3; i < 6; ++i){
	const string index = boost::lexical_cast<string>(i);
	TEST_ASSERT( curves.add(string("old")+index,(VectorXd)(Z1.row(i)),h,0,"o") );
	TEST_ASSERT( curves.add(string("new")+index,(VectorXd)(Z2.row(i)),h,0) );
  }
  TEST_ASSERT ( curves.write("tt2.py") );

  curves.clear();
  for (int i = 6; i < 10; ++i){
	const string index = boost::lexical_cast<string>(i);
	TEST_ASSERT( curves.add(string("old")+index,(VectorXd)(Z1.row(i)),h,0,"o") );
	TEST_ASSERT( curves.add(string("new")+index,(VectorXd)(Z2.row(i)),h,0) );
  }
  TEST_ASSERT ( curves.write("tt3.py") );
}

BOOST_AUTO_TEST_SUITE_END()
