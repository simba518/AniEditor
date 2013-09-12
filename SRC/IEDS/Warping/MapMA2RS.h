#ifndef _MAPMA2RS_H_
#define _MAPMA2RS_H_

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
using namespace Eigen;

namespace LSW_WARPING{
  
  /**
   * @class MapMA2RS Provide the matrices that map the Modal-Analysis
   * Coordinates to the Rotation-Strain Coordinates:
   * y = (PGW)z, where y are the RS coordinates and z are the MA coordinates, and
   * P: is a matrix that map the deformation gradient (GWz) to the RS coordinates.
   * G: the descrete deformation gradient operator, which can be obtained by
   *    using the class DefGradOperator.
   * W: the modal matrix, i.e., the eigen vectors.
   * @see the Study Record "Map Between Modal Coordinates and RS Coordinates".
   */
  class MapMA2RS{
	
  public:
	/**
	 * compute the map matrix PGW(e.g, P*G*W), by providing the descrete
	 * deformation operator G and Modal Matrix W.
	 * 
	 */
	static void computeMapMatPGW(const SparseMatrix<double> &G, const MatrixXd &W, MatrixXd &PGW);
	
	/**
	 * compute the map matrix P, where N is the number of the tetrahedron, and
	 * the result is returned in P.
	 * 
	 */
	static void computeMapMatP(const int N, SparseMatrix<double> &P);
	
  };
  
}//end of namespace

#endif /*_MAPMA2RS_H_*/
