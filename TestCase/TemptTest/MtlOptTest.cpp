#include <boost/test/unit_test.hpp>
#include <boost/filesystem.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <MtlOptEnergyAD.h>
#include "MtlOptModel.h"
using namespace Eigen;
using namespace UTILITY;
using namespace boost;
using namespace boost::unit_test::framework;

const int max_outter_it = 10000;
const double outer_tol = 1e-5;
const int save_per_it = 50;
string init_file="/home/simba/Workspace/AnimationEditor/Data/beam/mtlopt_cen_keyW.ini";

class MtlOptADTestSuite{
  
public:
  MtlOptADTestSuite(const string test_case_n):test_case_name(test_case_n){
	if(master_test_suite().argc >= 2){
	  init_file = master_test_suite().argv[1];
	}
	cout<< "init file: " << init_file << endl;
	cout<< "test opt method: " << test_case_n << endl;

	model= pMtlOptModel(new MtlOptModel(init_file));
	dataM = pMtlDataModel(new MtlDataModel);
	model->initMtlData(*dataM);
	model->print();
  }
  void addOptimizer(pMtlOptimizer opt){optimizer.push_back(opt);}
  void save(){

	filesystem::path inf_path = filesystem::path(init_file);
	const string dir = inf_path.parent_path().string();
	const string fname = inf_path.filename().string();
	const string f = dir+"/tempt/"+fname+test_case_name;

	TEST_ASSERT( EIGEN3EXT::write(f+"OptZ.b",dataM->Z) );
	TEST_ASSERT( EIGEN3EXT::write(f+"OptK.b",dataM->K) );
  }
  void printZ(){

	SelfAdjointEigenSolver<MatrixXd> eigenK(dataM->K);
	const MatrixXd &U = eigenK.eigenvectors();
	MatrixXd Z = dataM->Z;
	MatrixXd newZ = (U.transpose()*Z).transpose();
	cout << "curve_new_z: "<<newZ.rows()<<" "<<newZ.cols()<<" ";
	cout << Map<VectorXd >(&newZ(0,0),newZ.size()).transpose() << endl;
  }
  void solve(){

	for (int i = 0; i < max_outter_it; ++i){

	  const MatrixXd oldZ = dataM->Z;
	  for (int p = 0; p < optimizer.size(); ++p)
		optimizer[p]->optimize();

	  cout << "OUTTER ITERATION " << i << " --------" << endl;
	  dataM->print();
	  if(oldZ.size() == dataM->Z.size()){
		const double err = (oldZ-dataM->Z).norm()/dataM->Z.norm();
		cout << "OUTTER DIFFERRENCE: " << err << endl;
		if(err < outer_tol) break;
	  }

	  if (i>0 && i%save_per_it==0 ){
		printZ();
		save();
	  }
	}
	printZ();
	save();
  }
  pMtlDataModel getDataModel(){return dataM;}

private:
  const string test_case_name;
  vector<pMtlOptimizer> optimizer;
  pMtlOptModel model;
  pMtlDataModel dataM;
};

BOOST_AUTO_TEST_SUITE(MtlOptTest)

BOOST_AUTO_TEST_CASE(Opt_Z){

  MtlOptADTestSuite mtlOpt("Opt_Z");
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_ZAtA_simu){

  MtlOptADTestSuite mtlOpt("Opt_ZAtA_simu");
  mtlOpt.addOptimizer(pMtlOptimizer(new ZAtAOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_Z_Lambda){

  MtlOptADTestSuite mtlOpt("Opt_Z_Lambda");
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.addOptimizer(pMtlOptimizer(new LambdaOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_Z_AtA){

  MtlOptADTestSuite mtlOpt("Opt_Z_AtA");
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.addOptimizer(pMtlOptimizer(new AtAOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_Z_akam){

  MtlOptADTestSuite mtlOpt("Opt_Z_akam");
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.addOptimizer(pMtlOptimizer(new akamOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_Z_AtA_akam){

  MtlOptADTestSuite mtlOpt("Opt_Z_AtA_akam");
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.addOptimizer(pMtlOptimizer(new AtAakamOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_SUITE_END()
