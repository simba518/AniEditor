#include <iostream>
#include <vector>
#include <TetMesh.h>
#include <Mesquite_all_headers.hpp>
using namespace UTILITY;
using namespace Mesquite;
using namespace std;

int main(int argc, char *argv[]){

  TetMesh tetMesh;
  if(argc <= 1){
	cout << "usage: ./TetOpt input tet file (*.abq)" << endl;
	return -1;
  }

  const string tetFileName = argv[1];
  if(!tetMesh.load(tetFileName)){
	cout << "failed to load tetmesh from " << tetFileName << endl;
	return -2;
  }

  vector<double> nodes;
  vector<long unsigned int> tets;
  cout << tetMesh.tets().size() << endl;

  tetMesh.nodes(nodes);
  tetMesh.tets(tets);
  assert(nodes.size()>=4);
  assert(nodes.size()%3==0);
  assert(tets.size()>=1);
  assert(tets.size()%4==0);

  cout << tets.size()/4 << endl;

  vector<int> fixNode(nodes.size()/3, 0);
  vector<int> surface;
  tetMesh.surface(surface);
  for (size_t i = 0; i < surface.size(); ++i){
	assert(surface[i]>=0);
	assert(surface[i]<fixNode.size());
    fixNode[surface[i]] = 1;
  }
  
  MsqError err;
  ArrayMesh mesh(3,nodes.size()/3,&(nodes[0]),&fixNode[0],tets.size()/4,TETRAHEDRON,&tets[0]);

  ShapeImprover optimizer;
  // LaplaceWrapper optimizer;
  optimizer.run_instructions(&mesh, err);
  if (err) {
  	cerr << "ERROR: quality improve fail: "<< endl << err << endl;
  }

  tetMesh.reset(nodes,tets);
  if(!tetMesh.write(tetFileName+".abq")){
	cout << "ERROR: failed to save the abq file to "<< tetFileName+".abq" << endl;
  }
  if(!tetMesh.writeVTK(tetFileName+".vtk")){
	cout << "ERROR: failed to save the vtk file to "<< tetFileName+".vtk" << endl;
  }
  return 0;
}
