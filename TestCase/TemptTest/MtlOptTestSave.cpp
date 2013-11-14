#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MatrixIO.h>
#include "MtlOptModel.h"
using namespace Eigen;
using namespace UTILITY;
using namespace EIGEN3EXT;
using namespace boost::unit_test::framework;

BOOST_AUTO_TEST_SUITE(MtlOptTestSave)

BOOST_AUTO_TEST_CASE(saveAll){

  const string d = "/home/simba/Workspace/AnimationEditor/Data/beam/";
  string init_file = d + "mtlopt_cen_keyW.ini";
  string opt_z_str = d + "/tempt/mtlopt_cen_keyW.iniOpt_ZOptZ.b";
  string opt_k_str = d + "/tempt/mtlopt_cen_keyW.iniOpt_ZOptK.b";
  string method = "Opt_Z";

  if(master_test_suite().argc >= 5){
	init_file = master_test_suite().argv[1];
	opt_z_str = master_test_suite().argv[2];;
	opt_k_str = master_test_suite().argv[3];;
	method = master_test_suite().argv[4];;
  }
  
  MtlOptModel model(init_file);
  MatrixXd Z;
  MatrixXd K;
  TEST_ASSERT( load(opt_z_str,Z) );
  TEST_ASSERT( load(opt_k_str,K) );
  ASSERT_EQ ( K.rows(), K.cols() );
  ASSERT_EQ ( Z.rows(), K.cols() );
  ASSERT_GT ( K.size(),0 );
  SelfAdjointEigenSolver<MatrixXd> eigenK(K);
  const MatrixXd &U = eigenK.eigenvectors();
  
  boost::filesystem::path inf_path = boost::filesystem::path(init_file);
  const string dir = inf_path.parent_path().string();
  const string fname = inf_path.filename().string();
  const string mesh_dir_f = dir+"/tempt/mesh/"+fname+method;

  model.saveMesh(Z,mesh_dir_f+"output_mesh_");
  
  const MatrixXd newZ = U*Z;
  MatrixXd zi = newZ;
  for (int i = 0; i < newZ.rows() && i < 4; ++i){
    zi.setZero();
    zi.row(i) = newZ.row(i);
    model.saveMesh(zi,mesh_dir_f+"mode_"+TOSTR(i)+"_f");
  }
}

BOOST_AUTO_TEST_SUITE_END()
