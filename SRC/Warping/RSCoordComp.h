#ifndef _RSCOORDCOMP_H_
#define _RSCOORDCOMP_H_

#include <boost/shared_ptr.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <TetMesh.h>
using namespace Eigen;
using namespace UTILITY;

typedef std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > VectorV3;
typedef std::vector<Eigen::Matrix<double,6,1>,Eigen::aligned_allocator<Eigen::Matrix<double,6,1> > > VectorV6;

namespace LSW_WARPING{
  
  /**
   * @class RSCoordComp compute the RS coordinates from Euler coordinates.
   * 
   */
  class RSCoordComp{
	
  public:
	/** 
	 * Construct the RS coordinates with warping considered, where p is usually
	 * are linear displacements in Euler space with large artifacts, and the
	 * resulted RS coordinates represents the warped good-looking results. This
	 * function is used to compute the RS coordinates of the editing result "p"
	 * in the paper: "Interactive Editing of Deformable Simulations, Barbic
	 * sig2012".
	 * 
	 * @param G the deformation gradient operator calculated by DefGradOperator.
	 * @param p the displacements in Euler space.
	 * @param y the RS-coordinates.
	 */
	static void constructWithWarp(const SparseMatrix<double> &G, const VectorXd &p, VectorXd &y);

	/**
	 * Construct the RS coordinates without warping considered, where u is
	 * usually are the good-looking displacements in Euler space. And
	 * reconstruction from the result RS coordinates y will give results
	 * similarly to u. This function is used to compute the RS coordinates of
	 * the refrence sequence "\bar{u}" in the paper: "Interactive Editing of
	 * Deformable Simulations, Barbic sig2012".
	 * 
	 */
	static void constructWithoutWarp(const SparseMatrix<double> &G, const VectorXd &u, VectorXd &y);

	// given the RS coordinates y, compute the rotation vector w as the average
	// of the tetrahedron this nodes connected to.
	static void RS2NodeRotVec(pTetMesh_const tetmesh,const VectorXd &y, VectorV3 &w);

	// given the RS coordinates y, compute the sysmetric part of the RS
	// coordinates as the average of the tetrahedron this nodes connected to.
	static void RS2NodeRSVec(pTetMesh_const tetmesh,const VectorXd &y, VectorV6 &s);

  protected:
	/**
	 * compute the RS coordinates yj of one tetrahedron j with respect to the
	 * deformation gradient mj. Here, yj is a 9x1 vector, and mj is a 9x1 vector
	 * represents a row-major 3x3 matrix.
	 * 
	 */
	static void constructWithWarp_OneTet(const double *mj, double *yj);
	static void constructWithoutWarp_OneTet(const double *mj, double *yj);
	
	
  };
  
  typedef boost::shared_ptr<RSCoordComp> pRSCoordComp;
  
}//end of namespace

#endif /*_RSCOORDCOMP_H_*/
