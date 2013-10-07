#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>
#include <UnitTestAssert.h>
#include <DrawCurves.h>
#include <AuxTools.h>
#include <MapMA2RS.h>
#include <MeshVtkIO.h>
#include <TetMeshEmbeding.h>
#include <MASimulatorAD.h>
#include "ElasticMtlOpt.h"
#include "HarmonicOscillator.h"
#include "DFT.h"
using namespace Eigen;
using namespace LSW_SIM;
using namespace ANI_EDIT;
using namespace UTILITY;

BOOST_AUTO_TEST_SUITE(ElasticMtlTest)

BOOST_AUTO_TEST_CASE(computeTestUrs){
  
  ElasticMtlOpt mtlOpt(0.1f,30);
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  TEST_ASSERT(mtlOpt.loadVolMesh(data+"/sim-mesh.hinp"));
  TEST_ASSERT(mtlOpt.loadAniSeq(data+"/tempt/UstrCon10.b"));
  mtlOpt.compute();
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
  TetMeshEmbeding volobj;
  TEST_ASSERT( volobj.loadObjMesh(data+"beam.obj"));
  TEST_ASSERT( volobj.loadTetMesh(data+"sim-mesh.hinp"));
  TEST_ASSERT( volobj.loadWeights(data+"interp-weights.txt"));

  vector<int> fixed_nodes;
  TEST_ASSERT( loadVec(data+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
  
  RS2Euler rs2euler;
  rs2euler.setTetMesh(volobj.getTetMesh());
  // rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.fixBaryCenter();
  rs2euler.precompute();

  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(volobj.getTetMesh(),G) );
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
	cout<<"d: "<<(z-mtlopt._zrs.col(i)).norm()/mtlopt._zrs.col(i).norm()<<endl;
	Z.push_back(z);

  	TEST_ASSERT ( rs2euler.reconstruct(y,u) );
  	volobj.interpolate(u);
  	const string fname=data+"tempt/meshes/objnocon_"+TOSTR(i)+".vtk";
  	TEST_ASSERT (MeshVtkIO::write(volobj.getObjMesh(),fname));
	U.push_back(u);
  }
  
  const string Ustr = data+"/tempt/Urs"+TOSTR(r)+".b";
  TEST_ASSERT( write(Ustr,U));

  const string Zstr = data+"/tempt/Zrs"+TOSTR(r)+".b";
  TEST_ASSERT( write(Zstr,Z));

  const string Ystr = data+"/tempt/Yrs"+TOSTR(r)+".b";
  TEST_ASSERT( write(Ystr,Y));
  
  TEST_ASSERT( write(data+"/tempt/u0.b",U[0] ));
}

void updateRSCoords(vector<VectorXd> &Y){
  
  assert_gt(Y.size(),0);
  REP (e,Y[0].size()/9){
	for (int f = Y.size()-2; f >= 0; --f){
	  const Vector3d &w2 = Y[f+1].segment(e*9,3);
	  const Vector3d &w1 = Y[f].segment(e*9,3);
	  if(((w1.normalized().transpose())*w2.normalized())(0,0) < -0.9f ){
	  	const double theta = w1.norm();
	  	Y[f].segment(e*9,3) = (w1.normalized())*(theta-2.0f*M_PI);
	  }
	}
  }
}

BOOST_AUTO_TEST_CASE(checkRSCon){

  // load data
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  TetMeshEmbeding volobj;
  TEST_ASSERT( volobj.loadObjMesh(data+"beam.obj"));
  TEST_ASSERT( volobj.loadTetMesh(data+"sim-mesh.hinp"));
  TEST_ASSERT( volobj.loadWeights(data+"interp-weights.txt"));

  vector<int> fixed_nodes;
  TEST_ASSERT( loadVec(data+"/con_nodes.bou",fixed_nodes,UTILITY::TEXT) );
  
  RS2Euler rs2euler;
  rs2euler.setTetMesh(volobj.getTetMesh());
  rs2euler.setFixedNodes(fixed_nodes);
  rs2euler.precompute();

  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(volobj.getTetMesh(),G) );

  const int r = 5;
  const double h = 0.1f;
  MatrixXd tW,PGW;
  TEST_ASSERT ( load(data+"eigen_vectors80.b",tW) );
  const MatrixXd W = tW.leftCols(r);
  MapMA2RS::computeMapMatPGW(G,W,PGW);
  const MatrixXd invPGW = PseudoInverse(PGW);

  VectorXd u1,u2;
  vector<VectorXd> Y1,Y2;
  TEST_ASSERT( load(data+"/tempt/YrsCon10.b",Y1) );

  VectorXd eigenval;
  TEST_ASSERT ( load(data+"eigen_values80.b",eigenval) );
  cout<< "correct eigen values:\n" << eigenval.segment(0,r) << endl;

  // compute Y
  INFO_LOG("compute Y");
  const int T = Y1.size();
  MatrixXd Z1(r,T),Z2(r,T);
  REP (i,Y1.size()){
  	const VectorXd &y1 = Y1[i];
  	TEST_ASSERT ( rs2euler.reconstruct(y1,u1) );
  	VectorXd y2;
  	LSW_WARPING::RSCoordComp::constructWithoutWarp(G,u1,y2);
  	TEST_ASSERT ( rs2euler.reconstruct(y2,u2) );
  	Y2.push_back(y2);
  }

  updateRSCoords(Y2);
  REP(i,Y2.size()){
	Z1.col(i) = invPGW*Y1[i];
  	Z2.col(i) = invPGW*Y2[i];
	cout<<"i= "<<i<<",\t" << (Y1[i]-Y2[i]).norm() << endl;
  }

  // // analysis Y
  // INFO_LOG("analysis Y");
  // MatrixXd Ym(Y2[0].size(),Y2.size());
  // REP(f,Y2.size()) Ym.col(f) = Y2[f];
  // PythonFFT pfft;
  // PythonScriptDraw2DCurves<VectorXd> tetcurves;
  // // const int step = (Ym.rows()/9)/100;
  // const int step = 1;
  // REP(e,Ym.rows()/9){
  // 	if(e%step == 0){
  // 	  REP(i,9){
  // 		const int c = e*9+i;
  // 		const VectorXd curve = (VectorXd)Ym.row(c);
  // 		pfft.add(string("sig_")+TOSTR(e)+"_"+TOSTR(i),curve,h);
  // 		tetcurves.add(string("t_")+TOSTR(e)+"_"+TOSTR(i),curve,h);
  // 	  }
  // 	  // const string pname = string("./tempt2/fft")+TOSTR(e);
  // 	  // pfft.write(pname,string("./tempt2/figures/fft")+TOSTR(e));
  // 	  // const string cname = string("./tempt2/tetcurves")+TOSTR(e);
  // 	  // tetcurves.write(cname,string("./tempt2/figures/tetcurves")+TOSTR(e));
  // 	  // pfft.clear();
  // 	  // tetcurves.clear();
  // 	}
  // }
  // const string pname = string("./tempt2/fft");
  // pfft.write(pname,string("./tempt2/figures/fft")); 

  const int start = 20;

  // approximate
  INFO_LOG("approximate");
  ElasticMtlOpt mtlOpt(h,r);
  const MatrixXd subZ1 = Z1.rightCols(T);
  const MatrixXd subZ2 = Z2.rightCols(T);
  mtlOpt._zrs = subZ2;
  mtlOpt.computeK();
  mtlOpt._Wrs = PGW;
  mtlOpt.decomposeK();

  // // save curves
  // INFO_LOG("save curves");
  // PythonScriptDraw2DCurves<VectorXd> curves;
  // REP (i,Z1.rows()){
  // 	TEST_ASSERT(curves.add(string("o")+TOSTR(i+1),(VectorXd)(Z1.row(i)),h,0,"o"));
  // 	TEST_ASSERT(curves.add(string("n")+TOSTR(i+1),(VectorXd)(Z2.row(i)),h,0));
  // 	if((i+1)%3 == 0){
  // 	  TEST_ASSERT ( curves.write(string("t")+TOSTR((i+1)/3)+".py") );
  // 	  curves.clear();
  // 	}
  // }
  // curves.clear();

  // // save animations
  // INFO_LOG("save animations");
  // REP(mode,Z2.rows()-1){
  // 	const MatrixXd mZ = Z2.row(mode);
  // 	const MatrixXd mW = W.col(mode);
  // 	const MatrixXd minvW = invPGW.col(mode);
  // 	for (int i = 0; i < mZ.cols(); ++i){
  // 	  VectorXd y,u;
  // 	  const VectorXd p = mW*i*10.0f;
  // 	  RSCoordComp::constructWithWarp(G,p,y);
  // 	  TEST_ASSERT ( rs2euler.reconstruct(y,u) );
  // 	  TEST_ASSERT ( volobj.interpolate(&u[0]) );
  // 	  const string fname=data+"tempt/meshes/"+TOSTR(mode)+"obj_"+TOSTR(i)+".vtk";
  // 	  TEST_ASSERT (MeshVtkIO::write(volobj.getObjMesh(),fname));
  // 	}
  // }

  // // DFT
  // INFO_LOG("DFT");
  // cout << endl << endl << endl;
  // REP (i,Z1.rows()){
  // 	const VectorXd ci = Z1.row(i);
  // 	IOFormat OctaveFmt(10, 0, ", ", ";\n", "", "", "[", "]");
  // 	cout<<string("sig")+TOSTR(i)+"="<<ci.transpose().format(OctaveFmt)<<endl<<endl;
  // }
}

BOOST_AUTO_TEST_CASE(checkRSNoCon){

  // load data
  const string data = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  TetMeshEmbeding volobj;
  TEST_ASSERT( volobj.loadObjMesh(data+"beam.obj"));
  TEST_ASSERT( volobj.loadTetMesh(data+"sim-mesh.hinp"));
  TEST_ASSERT( volobj.loadWeights(data+"interp-weights.txt"));

  RS2Euler rs2euler;
  rs2euler.setTetMesh(volobj.getTetMesh());
  rs2euler.fixBaryCenter();
  rs2euler.precompute();

  SparseMatrix<double> G;
  TEST_ASSERT( LSW_WARPING::DefGradOperator::compute(volobj.getTetMesh(),G) );

  const int r = 10;
  const double h = 0.1f;
  MatrixXd tW,PGW;
  TEST_ASSERT ( load(data+"eigen_vectors_nocon80.b",tW) );
  const MatrixXd W = tW.leftCols(r);
  MapMA2RS::computeMapMatPGW(G,W,PGW);
  const MatrixXd invPGW = PseudoInverse(PGW);

  VectorXd u1,u2;
  vector<VectorXd> Y1,Y2;
  TEST_ASSERT( load(data+"/tempt/Yrs10.b",Y1) );

  VectorXd eigenval;
  TEST_ASSERT ( load(data+"eigen_values_nocon80.b",eigenval) );
  cout<< "correct eigen values:\n" << eigenval.segment(0,r).transpose() << endl;

  // compute Y
  INFO_LOG("compute Y");
  const int T = Y1.size();
  MatrixXd Z1(r,T),Z2(r,T);
  REP (i,Y1.size()){
  	const VectorXd &y1 = Y1[i];
  	TEST_ASSERT ( rs2euler.reconstruct(y1,u1) );
  	VectorXd y2;
  	LSW_WARPING::RSCoordComp::constructWithoutWarp(G,u1,y2);
  	TEST_ASSERT ( rs2euler.reconstruct(y2,u2) );
  	Y2.push_back(y2);
  }

  // updateRSCoords(Y2);
  REP(i,Y2.size()){
	Z1.col(i) = invPGW*Y1[i];
  	Z2.col(i) = invPGW*Y2[i];
	// cout << "f,\t" << (Y1[i]-Y2[i]).norm()/Y1[i].size() << endl;

	// int e = 0;
	// const VectorXd &s1 = Y1[i].segment(e*9,3);
	// const VectorXd &s2 = Y2[i].segment(e*9,3);
	// cout << s1.transpose() << endl;
	// cout << s2.transpose() << endl;
	
	// e = 100;
	// const VectorXd &s11 = Y1[i].segment(e*9,3);
	// const VectorXd &s21 = Y2[i].segment(e*9,3);
	// cout << s11.transpose() << endl;
	// cout << s21.transpose() << endl;
  }

  // // analysis Y
  // INFO_LOG("analysis Y");
  // MatrixXd Ym(Y2[0].size(),Y2.size());
  // REP(f,Y2.size()) Ym.col(f) = Y2[f];
  // PythonFFT pfft;
  // PythonScriptDraw2DCurves<VectorXd> tetcurves;
  // // const int step = (Ym.rows()/9)/100;
  // const int step = 1;
  // REP(e,Ym.rows()/9){
  // 	if(e%step == 0){
  // 	  REP(i,9){
  // 		const int c = e*9+i;
  // 		const VectorXd curve = (VectorXd)Ym.row(c);
  // 		pfft.add(string("sig_")+TOSTR(e)+"_"+TOSTR(i),curve,h);
  // 		tetcurves.add(string("t_")+TOSTR(e)+"_"+TOSTR(i),curve,h);
  // 	  }
  // 	  // const string pname = string("./tempt2/fft")+TOSTR(e);
  // 	  // pfft.write(pname,string("./tempt2/figures/fft")+TOSTR(e));
  // 	  // const string cname = string("./tempt2/tetcurves")+TOSTR(e);
  // 	  // tetcurves.write(cname,string("./tempt2/figures/tetcurves")+TOSTR(e));
  // 	  // pfft.clear();
  // 	  // tetcurves.clear();
  // 	}
  // }
  // const string pname = string("./tempt2/fft");
  // pfft.write(pname,string("./tempt2/figures/fft")); 

  const int start = 20;

  // approximate
  INFO_LOG("approximate");
  ElasticMtlOpt mtlOpt(h,r);
  const MatrixXd subZ1 = Z1.rightCols(T);
  const MatrixXd subZ2 = Z2.rightCols(T);
  MatrixXd realZ;
  TEST_ASSERT( load(data+"/tempt/Zrs10.b",realZ) );
  // Z1.row(6).normalize();
  // Z1.row(7).normalize();
  // Z1.row(8).normalize();
  // Z1.row(9).normalize();
  mtlOpt._zrs = Z2;
  mtlOpt.computeK();
  mtlOpt._Wrs = PGW;
  mtlOpt.decomposeK();

  // save curves
  INFO_LOG("save curves");
  PythonScriptDraw2DCurves<VectorXd> curves;
  REP (i,Z1.rows()){
  	TEST_ASSERT(curves.add(string("o")+TOSTR(i+1),(VectorXd)(Z1.row(i)),h,0,"o"));
  	TEST_ASSERT(curves.add(string("n")+TOSTR(i+1),(VectorXd)(Z2.row(i)),h,0));
  	if((i+1)%3 == 0){
  	  TEST_ASSERT ( curves.write(string("tempt/nocon_z")+TOSTR((i+1)/3)) );
  	  curves.clear();
  	}
  }
  curves.clear();

  // // save animations
  // INFO_LOG("save animations");
  // REPP(mode,0,9){
  // 	const MatrixXd mZ = Z1.row(mode);
  // 	const MatrixXd mW = W.col(mode);
  // 	const MatrixXd minvW = invPGW.col(mode);
  // 	for (int i = 0; i < mZ.cols(); ++i){
  // 	  VectorXd y,u;
  // 	  const VectorXd p = mW*mZ.col(i);
  // 	  RSCoordComp::constructWithWarp(G,p,y);
  // 	  TEST_ASSERT ( rs2euler.reconstruct(y,u) );
  // 	  TEST_ASSERT ( volobj.interpolate(&u[0]) );
  // 	  const string fname=data+"tempt/meshes/"+TOSTR(mode)+"noconobj_"+TOSTR(i);
  // 	  TEST_ASSERT (MeshVtkIO::write(volobj.getObjMesh(),fname+".vtk"));
  // 	}
  // }

  // // DFT
  // INFO_LOG("DFT");
  // cout << endl << endl << endl;
  // REP (i,Z1.rows()){
  // 	const VectorXd ci = Z1.row(i);
  // 	IOFormat OctaveFmt(10, 0, ", ", ";\n", "", "", "[", "]");
  // 	cout<<string("sig")+TOSTR(i)+"="<<ci.transpose().format(OctaveFmt)<<endl<<endl;
  // }
}

BOOST_AUTO_TEST_SUITE_END()
