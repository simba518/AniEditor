#include "StVKSimulator.h"
#include <ConMatrixTools.h>
#include <MatrixIO.h>
#include <Timer.h>
#include <CCD_DPF.h>
using namespace UTILITY;
using namespace SIMULATOR;

// load data from init file: mesh, damping, time step, constraints, forces
void StVKSimulator::loadInitFile(const string filename){
  
  if(!simulator->init(filename)){
	ERROR_LOG("failed to initialze the simulator.");
  }
  
  UTILITY::JsonFilePaser jsonf;
  if (jsonf.open(filename)){
	string fixed_nodes_str;
	if(jsonf.readFilePath("fixed_nodes",fixed_nodes_str)){
	  vector<int> fixed_nodes;
	  UTILITY::loadVec(fixed_nodes_str, fixed_nodes, UTILITY::TEXT);
	  cout << "number of fixed nodes: " << fixed_nodes.size() << endl;
	  setFixedNodes(fixed_nodes);
	}

	if(!jsonf.readFilePath("save_results_to", saveRlstTo,false)){
	  saveRlstTo = "tempt/simulated_results";
	}
	if(!jsonf.read("T",totalFrames)){
	  totalFrames = 200;
	}
	if(!jsonf.read("steps",steps)){
	  steps = 1;
	}

	vector<double> g;
	if(jsonf.read("gravity",g)){
	  assert_eq(g.size(),3);
	  gravity[0] = g[0];
	  gravity[1] = g[1];
	  gravity[2] = g[2];
	}else{
	  clearGravity();
	}

	jsonf.read("gravity_start", gravity_start);
	jsonf.read("gravity_stop", gravity_stop);
	jsonf.readVecFile("u0",u0);

  }else{
	ERROR_LOG("failed to open the initfile: " << filename);
  }
}

void StVKSimulator::initCollision(){

  CCD_DPF::getInstance()->addObject();
  CCD_DPF::getInstance()->prepare();
}

void StVKSimulator::simulate(){

  bool succ = true;
  recorded_U.clear();
  recorded_U.reserve(totalFrames);
  cout << "total steps: " << totalFrames << endl;
  cout << "gravity start: " << gravity_start << endl;
  cout << "gravity stop: " << gravity_stop << endl;

  Timer timer;
  timer.start();
  if (u0.size() > 0){
	assert_eq(u0.size(), stvkModel->dimension());
	simulator->setU0(u0);
  }

  for (int i = 0; i < totalFrames; ++i){

	cout << "step = " << i << endl;
	if (i >= gravity_start && i <= gravity_stop){
	  simulator->setExtForceForAllNodes(gravity[0], gravity[1], gravity[2]);
	}else{
	  double g[3] = {0,0,0};
	  simulator->setExtForceForAllNodes(g[0],g[1],g[2]);
	}

	CCD_DPF::getInstance()->update();
	CCD_DPF::checkCollision(,);
	CCD_DPF::computeForce();
	CCD_DPF::addCollisionForce(stvkModel->getTetMesh().get(),f_ext);

	bool succ = false;
	for (int i = 0; i < steps; ++i){
	  succ = simulator->forward();
	}

	recorded_U.push_back(simulator->getU());
	ERROR_LOG_COND("simulation failed, step i = " << i, succ);
  }

  timer.stop("total simulation time is: ");
}

void StVKSimulator::save(){
  
  if(!EIGEN3EXT::write(saveRlstTo+".b",recorded_U)){
	ERROR_LOG("failed to save the results U to " << saveRlstTo<<".b");
  }else{
	cout << "success to save the results U to: " << saveRlstTo<<".b" << endl;
  }

  const pTetMesh_const tetmesh = stvkModel->getTetMesh();
  assert(tetmesh);
  if(!tetmesh->writeVTK(saveRlstTo, recorded_U)){
	ERROR_LOG("failed to save the results mesh to " << saveRlstTo<<"*.vtk");
  }else{
	cout << "success to save the results mesh to: " << saveRlstTo<<"*.vtk" << endl;
  }
  
}

void StVKSimulator::setFixedNodes(const vector<int> &fixednodes){

  VectorXd uc(fixednodes.size()*3);
  uc.setZero();
  simulator->setUc(uc);
  
  VecT trip_C;
  const int n = stvkModel->dimension()/3;
  UTILITY::computeConM(fixednodes, trip_C, n);
  simulator->setConM(trip_C);
}
