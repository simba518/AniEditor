#ifndef _RS2EULER_H_
#define _RS2EULER_H_

#include <vector>
#include <set>
#include <boost/shared_ptr.hpp>
#include <TetMesh.h>
#include <Eigen/UmfPackSupport>
#include <eigen3/Eigen/Sparse>
using namespace std;
using namespace Eigen;
using namespace UTILITY;

namespace LSW_WARPING{

  /**
   * @class RS2Euler construct the Euler coordinates from the given
   * Rotation-Strain coordinates, and providing all of the functions to prepare
   * the construction, e.g. assemble and factor the left side coefficients
   * matrix, right side vectors.
   * 
   */
  class RS2Euler{
	
  public:
	void setTetMesh(pTetMesh_const tetmesh);
	void fixBaryCenter();
	void setFixedNodes(const vector<int> &fixed_nodes);
	void setFixedNodes(const set<int> &fixed_nodes);
	void setConNodes(const vector<set<int> > &c_nodes,const VectorXd &bcen_uc);
	void setConNodes(const vector<int> &c_nodes,const VectorXd &uc);
	bool precompute();	
	bool reconstruct(const VectorXd &y, VectorXd &u);

	const SparseMatrix<double> &get_G()const{
	  return G;
	}
	const SparseMatrix<double> &get_A()const{
	  return A;
	}
	const SparseMatrix<double> &get_L()const{
	  return L;
	}
	const SparseMatrix<double> &get_VG_t()const{
	  return VG_t;
	}
	const vector<double> &get_Sqrt_V()const{
	  return Sqrt_V;
	}
	pTetMesh_const getTetMesh()const{
	  return tetmesh;
	}
	void assemble_b(const VectorXd &y,VectorXd &b)const;


  protected:
	inline void constructVolumeMatrix(pTetMesh_const tetmesh, SparseMatrix<double> &V,vector<double> &Sqrt_V)const;
	inline void assembleA();
	inline int numFixedDofs()const{
	  return F.rows()+C.rows();
	}
	
  private:
	pTetMesh_const tetmesh; // tetrahedron mesh.

	/*
	 * A = | L C^t F^t|
	 *     | C  O  O  |
	 *     | F  O  O  |
	 */	
	SparseMatrix<double> A;
	SparseMatrix<double> L; // VG_t*VG
	SparseMatrix<double> F; // constraint matrix of fixed nodes
	SparseMatrix<double> C; // constraint matrix of constrained nodes
	VectorXd barycenter_uc;// barycenters of each group of the constrained nodes

	SparseMatrix<double> G; // deformation gradient operator
	SparseMatrix<double> VG_t; // (VG)^t
	vector<double> Sqrt_V; /// sqrt(V_i), where V_i is the volume of the i-th
						   /// tetrahedron.
	UmfPackLU<SparseMatrix<double> > A_solver;

  };
  
  typedef boost::shared_ptr<RS2Euler> pRS2Euler;
  
}//end of namespace

#endif /*_RS2EULER_H_*/
