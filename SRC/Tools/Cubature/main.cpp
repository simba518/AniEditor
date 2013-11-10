#include <iostream>
#include <eigen3/Eigen/Dense>
#include <LogConfig.h>
#include <JsonIniFReader.h>
#include <MatrixIO.h>
#include <volumetricMeshLoader.h>
#include "PerNodeCubature.h"
#include "ReducedRSCubature.h"
using namespace std;
using namespace Eigen;
using namespace LSW_UTILITY;

MatrixXd W,B;
VectorXd eigenvalues;
pVolumetricMesh tetmesh;
MatrixXd training_z;
MatrixXd training_u; // q(z) = B^t*u(z)
TrainingSet trainingZ;
VECTOR trainingQ;
double rel_err_tol;
double max_num_points;
double num_cands_per_iter;
double iters_per_fullnnls;
double num_samples_per_subtrain;
string save_selpoints;
string save_weights;
JsonIniFReader json_f;

bool loadInitData(int argc, char *argv[]){
  
  LogConfig::getInstance()->config();
  if(argc < 2){
	cout << "usage: ./Cubature init_filename" <<endl;
	return false;
  }
  const string ini_filename = argv[1];
  if (json_f.open(ini_filename)){

	json_f.read("rel_err_tol",rel_err_tol);
	json_f.read("max_num_points",max_num_points);
	json_f.read("iters_per_fullnnls",iters_per_fullnnls);
	json_f.read("num_samples_per_subtrain",num_samples_per_subtrain);
	json_f.read("num_cands_per_iter",num_cands_per_iter);
	json_f.readFilePath("save_selpoints",save_selpoints,false);
	json_f.readFilePath("save_weights",save_weights,false);

	json_f.read("eigen_vectors",W);
	json_f.read("eigen_values",eigenvalues);
	json_f.read("nonlinear_basis",B);
	string tetfile;
	if(json_f.readFilePath("vol_filename",tetfile)){
	  tetmesh = pVolumetricMesh(VolumetricMeshLoader::load(tetfile.c_str()));
	}else{
	  return false;
	}
  }
  return true;
}

void initTrainingSet(){
  
  if(!json_f.read("training_z",training_z)){

	// gen training set for each mode
	const int trainingPerMode = 20;
	const int numMode = W.cols();
	const int trainingMode = 20;
	training_z.resize(numMode,trainingPerMode*trainingMode);
	double scale = 100.0f;
	if(!json_f.read("mode_scale",scale)){
	  cout << "WARNNING: failed to read mode_scale, use 100." << endl;
	  scale = 100.0f;
	}

	for (int mode = 0; mode < trainingMode; ++mode){
	  assert_gt(eigenvalues[mode],0);
	  for (int t = 0; t < trainingPerMode; ++t){
		const double z_mode = (scale*t/trainingPerMode)/(sqrt(eigenvalues[mode]));
		training_z.col(mode*trainingPerMode+t).setZero();
		training_z.col(mode*trainingPerMode+t)[mode] = z_mode;
	  }
	}
  }
}

void runAndSave(WarpingCubature *pcubature){
  
  // compute 
  pcubature->run(trainingZ, 
				 trainingQ, 
				 rel_err_tol, 
				 max_num_points,
				 num_cands_per_iter,
				 iters_per_fullnnls,
				 num_samples_per_subtrain);

  // save results
  LSW_UTILITY::save(pcubature->getSelPoints(), save_selpoints,LSW_UTILITY::TEXT);
  LSW_UTILITY::save(pcubature->getWeights(), save_weights);

  // save as vtk format
  pcubature->saveAsVTK(save_selpoints+string(".vtk"));
}

int RunReducedRSCubature(){
  
  ReducedRSCubature cubature(W,B,tetmesh);
  initTrainingSet();

  if(!json_f.read("training_u",training_u)){
	cubature.initTrainingData(training_z,trainingZ, trainingQ);
  }else{
	const MatrixXd training_q = B.transpose()*training_u;
	cubature.convertTrainingData(training_z,training_q,trainingZ, trainingQ);
  }

  runAndSave(&cubature);

  // check the difference
  vector<VectorXd> z_for_check;
  string fixed_nodes_str;
  vector<int> fixed_nodes;
  const double vol = tetmesh->getTotalVolume();
  if (json_f.read("z_for_check",z_for_check) && 
	  json_f.readFilePath("con_nodes",fixed_nodes_str) &&
	  LSW_UTILITY::read(fixed_nodes,fixed_nodes_str,LSW_UTILITY::TEXT) ){

	cubature.setFixedNodes(fixed_nodes);
	for (size_t i = 0; i < z_for_check.size(); ++i){

	  const VectorXd &z = z_for_check[i];
	  const VectorXd rs_u = cubature.computeRSDisp(z);
	  const VectorXd u = B*(cubature.computeCubDisp(z));
	  cout << "RS diff REDUCED-RS: " << (u-rs_u).norm()/vol << endl;
	}
  }
  return 1;
}

int RunPerNodeCubature(){

  PerNodeCubature cubature(W,B,tetmesh);
  cubature.setNodeToCubature(10);
  initTrainingSet();
  cubature.initTrainingData(training_z,trainingZ, trainingQ);
  runAndSave(&cubature);
  return 1;
}

/**
 * given an init file including eigen vectors, eigen values, nonlinear basis, as
 * well as many other parameters, then it will calculate the cubatured
 * tetrahedrons and corresponding weights, and save them to files.
 */
int main(int argc, char *argv[]){

  if(!loadInitData(argc,argv)){
	cout << "failed to load init data.\n";
	return -1;
  }
  return RunReducedRSCubature();
}
