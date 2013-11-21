#ifndef _VOLAPPROXTRI_H_
#define _VOLAPPROXTRI_H_

#include <ElasticForceTetFullStVK.h>
#include <TetMeshEmbeding.h>
using namespace UTILITY;

namespace UTILITY{
  
  /**
   * @class VolApproxTri approximate the triangle mesh using volumetric mesh.
   * 
   * input:
   * 1. rest shape of tetrahedron mesh: r.
   * 2. rest shape of obj mesh: x0.
   * 3. keyframe of obj mesh: xk.
   * 4. interpolation weights.
   * 
   * output:
   * keyframe of tetrahedron mesh.
   * 
   **/
  class VolApproxTri{
	
  public:
	VolApproxTri(){
	  _tetMeshRest = pTetMesh(new TetMesh);
	  _objMeshRest = pObjmesh(new Objmesh);
	  _objMeshKey = pObjmesh(new Objmesh);
	  _penalty = 1.0f;
	}
	bool loadRestVolMesh(const string filename);
	bool loadRestObjMesh(const string filename);
	bool loadKeyObjMesh(const string filename);
	bool loadInterpWeights(const string filename);
	bool generateKeyVolMesh();
	bool saveAll(const string filename);

	void funGrad(const VectorXd &u,double &fun,VectorXd &grad);

  protected:
	void prepare();
	
  private:
	pTetMesh _tetMeshRest;
	pObjmesh _objMeshRest;
	pObjmesh _objMeshKey;
	TetMeshEmbeding _embed;
	VectorXd _volKeyU;
	double _penalty;

	VectorXd _volRestU;
	VectorXd _objKeyU;
	SparseMatrix<double> _A;
	pElasticForceTetFullStVK _elasticModel;
  };
  
}//end of namespace

#endif /*_VOLAPPROXTRI_H_*/
