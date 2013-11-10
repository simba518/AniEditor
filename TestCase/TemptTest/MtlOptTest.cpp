#include <boost/test/unit_test.hpp>
#include <UnitTestAssert.h>
#include <eigen3/Eigen/Dense>
#include <DrawCurves.h>
#include <MtlOptEnergyAD.h>
#include "MtlOptModel.h"
using namespace Eigen;
using namespace UTILITY;

class MtlOptADTestSuite{
  
public:
  MtlOptADTestSuite(){
	dataDir = "/home/simba/Workspace/AnimationEditor/Data/";
	maxOurterIt = 30;
	dataM = pMtlDataModel(new MtlDataModel);
  }
  MtlOptADTestSuite(const string d,
					const string initf, 
					const string name,
					const int outIt = 30,
					const string saveInputM = ""
					){

	dataDir = "/home/simba/Workspace/AnimationEditor/Data/";
	data = dataDir+d+"/";
	initFile = data+initf;
	maxOurterIt = outIt;
	
	outputMesh = data+"/tempt/mesh/"+initf+"_"+"_"+name+"_outPutMesh";
	curveZName = data+"/tempt/mesh/"+initf+"_"+name+"_curveZ";
	savePartialCon = data+"/tempt/mesh/"+initf+"_"+name+"_parcon";
	saveKeyframes = data+"/tempt/mesh/"+initf+"_"+name+"_keyf";
	saveModes = data+"/tempt/mesh/"+initf+"_"+name+"_modes";

	if (saveInputM.size() > 0) inputMesh = data+saveInputM;

	model= pMtlOptModel(new MtlOptModel(initFile));
	model->produceSimRlst(false);

	dataM = pMtlDataModel(new MtlDataModel);
	model->initMtlData(*dataM);
  }
  void addOptimizer(pMtlOptimizer opt){
	optimizer.push_back(opt);
  }
  void solve(){	

	for (int i = 0; i < maxOurterIt; ++i){

	  const MatrixXd oldZ = dataM->Z;
	  for (int p = 0; p < optimizer.size(); ++p){
		optimizer[p]->optimize();;
	  }

	  cout << "ITERATION " << i << "------------------------------------" << endl;
	  dataM->print();
	  if(oldZ.size() == dataM->Z.size()){
		const double err = (oldZ-dataM->Z).norm()/dataM->Z.norm();
		cout << "err: " << err << endl;
		if(err < 1e-3) break;
	  }
	}

	if (outputMesh.size() > 0) model->saveMesh(dataM->Z,outputMesh);
	if (inputMesh.size() > 0) model->saveMesh(model->Z,inputMesh);
	if (saveKeyframes.size() > 0 && model->Kz.size()>0) 
	  model->saveMesh(model->Kz,saveKeyframes);
	if (savePartialCon.size() > 0) model->saveUc(savePartialCon);

	PythonScriptDraw2DCurves<VectorXd> curves;
	const MatrixXd newZ = dataM->getRotZ();
	const MatrixXd newKZ = dataM->getUt()*model->Kz;
	for (int i = 0; i < model->Z.rows(); ++i){
	  const VectorXd &z2 = newZ.row(i).transpose();
	  curves.add(string("mode ")+TOSTR(i),z2,1.0f,0,"--");
	  const VectorXd &kz = newKZ.row(i).transpose();
	  VectorXd kid(model->Kid.size());
	  for (int j = 0; j < kid.size(); ++j)
		kid[j] = model->Kid[j];
	  curves.add(string("keyframes of mode ")+TOSTR(i),kz,kid,"o");
	}
	if (curveZName.size() > 0)
	  TEST_ASSERT( curves.write(curveZName) );

	const MatrixXd newW = model->W*(dataM->getUt().transpose());
	model->W = newW;
	MatrixXd zi = newZ;
	for (int i = 0; i < newZ.rows() && i < 4; ++i){
	  zi.setZero();
	  zi.row(i) = newZ.row(i);
	  model->saveMesh(zi,saveModes+TOSTR(i)+"_f");
	}
  }
  pMtlDataModel getDataModel(){
	return dataM;
  }

private:
  string dataDir;
  string data;
  string initFile;
  int maxOurterIt;
  string outputMesh;
  string curveZName;
  string inputMesh;
  string savePartialCon;
  string saveKeyframes;
  string saveModes;
  vector<pMtlOptimizer> optimizer;

  pMtlOptModel model;
  pMtlDataModel dataM;
};

const string model_name = "beam";
const string inf = "mtlopt_cen_keyf_swing_one_end.ini";

BOOST_AUTO_TEST_SUITE(MtlOptTest)

BOOST_AUTO_TEST_CASE(Opt_Z){

  MtlOptADTestSuite mtlOpt(model_name,inf,"Opt_Z",1);
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_Z_Lambda){

  MtlOptADTestSuite mtlOpt(model_name,inf,"Opt_Z_Lambda",500);
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.addOptimizer(pMtlOptimizer(new LambdaOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_Z_AtA){

  MtlOptADTestSuite mtlOpt(model_name,inf,"Opt_Z_AtA",50);
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.addOptimizer(pMtlOptimizer(new AtAOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_Z_akam){

  MtlOptADTestSuite mtlOpt(model_name,inf,"Opt_Z_akam",50);
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.addOptimizer(pMtlOptimizer(new akamOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_CASE(Opt_Z_AtA_akam){

  MtlOptADTestSuite mtlOpt(model_name,inf,"Opt_Z_AtA_akam",50);
  mtlOpt.addOptimizer(pMtlOptimizer(new ZOptimizer(mtlOpt.getDataModel())));
  mtlOpt.addOptimizer(pMtlOptimizer(new AtAakamOptimizer(mtlOpt.getDataModel())));
  mtlOpt.solve();
}

BOOST_AUTO_TEST_SUITE_END()
