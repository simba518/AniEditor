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
	  _penalty = 0.000001f;
	}
	void setEmbedingMesh(pTetMeshEmbeding embed){
	  _embed = embed;
	}
	void setElasticPenalty(const double p){
	  assert_gt(p,0.0f);
	  _penalty = p;
	}
	bool generateKeyVolMesh(pObjmesh_const objMeshKey);
	void funGrad(const VectorXd &u,double &fun,VectorXd &grad);
	void fun(const VectorXd &u,double &fun);
	const VectorXd &getVolKeyU()const{
	  return _volKeyU;
	}

  protected:
	void prepare(pObjmesh_const objMeshKey);
	bool alglibErrorReport(const int code)const;
	
  private:
	pTetMeshEmbeding _embed;
	VectorXd _volKeyU;
	double _penalty;
	VectorXd _volRestU;
	VectorXd _objKeyU;
	SparseMatrix<double> _A;
	pElasticForceTetFullStVK _elasticModel;
  };
  
}//end of namespace

#endif /*_VOLAPPROXTRI_H_*/
