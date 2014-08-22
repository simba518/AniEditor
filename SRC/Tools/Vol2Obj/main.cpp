#include <TetMeshEmbeding.h>
using namespace UTILITY;

int main(int argc, char *argv[]){

  if (argc < 3){
	cout << "error: incorrect number of arguments\n";
	cout << "ussage: " << "vol2obj tetfile objfile [saveto]" << endl;
  }

  const string tetfile = argv[1];
  const string objfile = argv[2];
  const string saveto = argc>3 ? argv[3]:tetfile+".interp";

  TetMeshEmbeding embeding;
  if( !embeding.loadTetMesh(tetfile) || !embeding.loadObjMesh(objfile)){
	return -1;
  }

  embeding.buildInterpWeights();
  if (embeding.writeWeights(saveto)){
	cout << "success to save the interpolate weights to " << saveto << endl;
  }else{
	cout << "error: failed to save the interpolate weights to " << saveto << endl;
  }
  return 0;
}
