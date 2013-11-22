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
	  _penalty = 0.000001f;
	}
	bool loadRestVolMesh(const string filename);
	bool loadRestObjMesh(const string filename);
	bool loadKeyObjMesh(const string filename);
	bool loadInterpWeights(const string filename);
	void setElasticPenalty(const double p){
	  assert_gt(p,0.0f);
	  _penalty = p;
	}
	bool generateKeyVolMesh();
	bool generateKeyVolMeshNumDiff();
	bool saveAll(const string filename);

	double penalty()const{
	  return _penalty;
	}
	void funGrad(const VectorXd &u,double &fun,VectorXd &grad);
	void fun(const VectorXd &u,double &fun);

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
