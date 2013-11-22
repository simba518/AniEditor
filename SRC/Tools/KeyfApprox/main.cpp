#include "VolApproxTri.h"
using namespace UTILITY;

int main(int argc, char *argv[]){
    
  if(argc <= 5){
	cout << "usage: KeyfApprox rest_tet rest_obj key_obj interp_weights save_to [elastic_penalty]" << endl;
	exit(0);
  }
  
  VolApproxTri approx;

  while(true){

	// load mesh and interp weights
	if(!approx.loadRestVolMesh(argv[1])) break;
	if(!approx.loadRestObjMesh(argv[2])) break;
	if(!approx.loadKeyObjMesh(argv[3])) break;
	if(!approx.loadInterpWeights(argv[4])) break;
	if(argc > 6){
	  approx.setElasticPenalty(atof(argv[6]));
	}else{
	  cout << "defualt elastic penalty will be used: "<< approx.penalty() << endl;
	}
	cout << "success to load all data, approximation is begining....."<<endl;

	// approximation
	if(!approx.generateKeyVolMeshNumDiff()){
	  cout << "error: failed to generate the keyframe of the volumetric mesh." << endl;
	}else{
	  cout << "success to approximate the obj keyframes with a tetrahedron mesh."<<endl;
	}

	// save
	if(approx.saveAll(argv[5])){
	  cout << "success to save all data to: " << argv[5] << endl;
	}else{
	  cout << "error: failed to save the data." << endl;
	}
	
	break;
  }

  return 0;
}
