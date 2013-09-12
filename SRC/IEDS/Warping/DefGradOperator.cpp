#include <stdlib.h>
#include <Log.h>
#include "DefGradOperator.h"
using namespace LSW_WARPING;

bool DefGradOperator::compute(pVolumetricMesh_const tetmesh,SparseMatrix<double> &G){

  assert(tetmesh != NULL);
  bool succ = true;
  if (tetmesh->numElementVertices() != 4){
	ERROR_LOG("only tetrahedron mesh supported.");
	succ = false;
  }else{

	const int elem_num = tetmesh->numElements();
	const int node_num = tetmesh->numVertices();
	const int nonzeros = 9*12*elem_num;

	vector<T> G_triplets;
	G_triplets.reserve(nonzeros);

	double V[4][3];
	double Ge[9][12];
	memset(&V[0][0],0,4*3*sizeof(double));
	memset(&Ge[0][0],0,9*12*sizeof(double));

	for (int e = 0; e < elem_num; ++e){

	  tetmesh->vertexOfTetEle(e,V);
	  ComputeGeAutoDiff::getInstance()->compute(V,Ge);
	  const int *elem_v = tetmesh->element(e)->vertices();
	  assembleGe2G(elem_v,e,Ge,G_triplets);
	}

	G.resize( elem_num*9, node_num*3 );
	G.reserve( nonzeros );
	G.setFromTriplets( G_triplets.begin(),G_triplets.end() );
	G.makeCompressed();
  }
  return succ;
}
