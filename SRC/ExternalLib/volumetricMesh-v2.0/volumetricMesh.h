/*************************************************************************
 *                                                                       *
 * "volumetricMesh" library , Copyright (C) 2007 CMU, 2009 MIT           *
 * All rights reserved.                                                  *
 *                                                                       *
 * Code author: Jernej Barbic                                            *
 * http://www.jernejbarbic.com/code                                      *
 * Research: Jernej Barbic, Doug L. James, Jovan Popovic                 *
 * Funding: NSF, Link Foundation, Singapore-MIT GAMBIT Game Lab          *
 * Version 2.0                                                           *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of the BSD-style license that is            *
 * included with this library in the file LICENSE.txt                    *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the file     *
 * LICENSE.TXT for more details.                                         *
 *                                                                       *
 *************************************************************************/

#ifndef _VOLUMETRICMESH_H_
#define _VOLUMETRICMESH_H_

/*
  This abstract class can store a generic volumetric 3D mesh. 
  It stores the mesh geometric information (elements and vertices), 
  and also the material parameters of each individual mesh element 
  (Young's modulus, Poisson ratio, mass density). This is done by 
  organizing elements with the same material parameters into a "solid section".
  The class supports several geometric queries and interpolation to 
  an embedded triangle mesh ("Free-Form Deformation").  It also 
  supports exporting the mesh to an .ele or .node format (the format 
  used, e.g., by the TetGen mesh generation package).  Derived classes 
  are the tetmesh (general tetrahedral meshes), and sqmesh 
  (axis-aligned voxel meshes). See description in sqmesh.h and tetmesh.h.

  All quantities are 0-indexed.

  See volumetricMeshInfo.cpp for a simple driver example.
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <vector>
#include <set>
#include <map>
#include <boost/shared_ptr.hpp>
#include <assertext.h>
#include "minivector.h"
#include "meshElement.h"
#include "elementSet.h"
#include "solidSection.h"
#include "material.h"

class VolumetricMesh
{
public:

  // constructs a mesh from the given vertices and elements, with a single solid
  // section and material "vertices" is double-precision array of length
  // 3xnumVertices .  "elements" is an integer array of length numElements x
  // numElementVertices
  VolumetricMesh(int numVertices, double * vertices,int numElements, int numElementVertices, 
				 int * elements,double E, double nu, double density); 
  
  // creates a mesh consisting of the specified element subset of the given
  // sqmesh if vertexMap is non-null, it also returns a renaming datastructure:
  // vertexMap[big mesh vertex] is the vertex index in the subset mesh
  VolumetricMesh(VolumetricMesh &mesh, int numElements,
				 int * elements, std::map<int,int> * vertexMap = NULL);

  virtual ~VolumetricMesh();

  // looks into the file and returns the element type (a string) of the
  // volumetric mesh in the file returns a NULL pointer if no type information
  // found the returned pointer to the string is allocated inside the routine
  // (unless NULL pointer is returned)
  static char * getElementType(const char * filename, char ** elementType); 

  // calls the derived class to identify itself
  virtual const char * getElementType() = 0; 

  // saves the mesh to a text file (a format similar to Abaqus is used; only
  // geometry is saved; not material information)
  virtual int save(char * filename) = 0; 

  // exports the elements to an .ele and .node file if includeSolidSections=1,
  // an extra column is added, identifying the solid section of each element
  int exportToEle(char * baseFilename, int includeSolidSections=0); 

  // vertex access
  inline int numVertices() const{ return numVertices_; }
  inline const Vec3d * vertex(int i) const{ return vertices_[i]; }//added by lsw
  inline void vertex(int i,double v[3]) const{ //added by lsw
	assert(i>=0 && i< this->numVertices());
	const Vec3d &v3 = *(this->vertex(i));
	v[0] = v3[0];
	v[1] = v3[1];
	v[2] = v3[2];
  }
  void vertices(std::vector<double> &v)const; // added by lsw

  inline Vec3d * vertex(int i) { return vertices_[i]; }
  inline Vec3d * vertex(int element, int vertex) { return vertices_[elements_[element]->vertex(vertex)]; }
  inline const Vec3d * vertex(int element, int vertex) const { return vertices_[elements_[element]->vertex(vertex)]; }
  inline int vertexIndex(int element, int vertex)const { return elements_[element]->vertex(vertex); }
  inline Vec3d ** vertices() { return vertices_;} // advanced, internal datastructure

  // element access
  inline int numElements() const{ return numElements_; }
  Vec3d elementCenter(int el)const;
  inline MeshElement * element(int i) { return elements_[i]; }
  inline const MeshElement * element(int i)const { return elements_[i]; } // add by lsw 
  inline MeshElement ** elements() { return elements_; }
  inline int numElementVertices() const{ return numElementVertices_; }
  inline const Vec3d &vertexOfElement(int el,int i)const{
	assert_ge(i,0);
	assert_ge(el,0);
	assert_gt(numElementVertices(),i);
	assert_gt(numElements(),el);
	const Vec3d *v = vertex(element(el)->vertex(i));
	assert (v!=NULL);
	return *v;
  } // add by lsw
  inline void vertexOfTetEle(int el,double V[4][3])const{

	const Vec3d & v0 = vertexOfElement(el,0);
	const Vec3d & v1 = vertexOfElement(el,1);
	const Vec3d & v2 = vertexOfElement(el,2);
	const Vec3d & v3 = vertexOfElement(el,3);
	V[0][0] = v0[0];   V[0][1] = v0[1];   V[0][2] = v0[2];  
	V[1][0] = v1[0];   V[1][1] = v1[1];   V[1][2] = v1[2];  
	V[2][0] = v2[0];   V[2][1] = v2[1];   V[2][2] = v2[2];  
	V[3][0] = v3[0];   V[3][1] = v3[1];   V[3][2] = v3[2];  
  } // add by lsw, only for tetrahedral mesh.
  
  /// return the barycenters of a set of groups of vertices. the result is
  /// stored in bary_p, a 3x1 vector for each group.
  /// @note bary_p should be pre-allocated, and len(bary_p) >= g.size()*3.
  void barycentersOfGroup(const std::vector<std::set<int> >&g, double *bary_p)const;// add by lsw

  /// return the barycenters of a groups of vertices. 
  void barycentersOfGroup(const std::vector<int>&one_group, double bary_p[3])const;// add by lsw

  /// return the barycenters of a groups of vertices. 
  void barycentersOfGroup(const std::set<int>&one_group, double bary_p[3])const;// add by lsw

  // materials
  inline int numMaterials() { return numMaterials_; }
  inline int numElementSets() { return numElementSets_; }
  inline int numSolidSections() { return numSolidSections_; }
  inline Material * material(int i) { return materials_[i]; }
  inline ElementSet * elementSet(int i) { return elementSets_[i]; }
  inline SolidSection * solidSection(int i) { return solidSections_[i]; }
  // sets i-th material to the given material
  void setMaterial(int i, Material & material);

  // erases all materials and creates a single material for the entire mesh
  void setSingleMaterial(double E, double nu, double density);

  // convert this mesh to a mesh containing the specified elements (i.e., delete
  // the mesh elements not on the given list of elements)
  void subsetMesh(std::set<int> & elements, int removeIsolatedVertices=1);

  // mesh 1-neighborhood queries
  void getVerticesInElements(std::vector<int> & elements, std::vector<int> & vertices);
  void getElementsTouchingVertices(std::vector<int> & vertices, std::vector<int> & elements);
  void getVertexNeighborhood(std::vector<int> & vertices, std::vector<int> & neighborhood);
 
  // center of mass and inertia tensor
  virtual double getElementVolume(int el)const = 0;
  double getTotalVolume()const{
	double vol = 0.0f;
	for (int i = 0; i < numElements(); ++i){
	  vol += getElementVolume(i);
	}
	return vol;
  }

  // returns the inertia tensor of a single element, around its center of mass,
  // with unit density
  virtual void getElementInertiaTensor(int el, Mat3d & inertiaTensor) = 0; 

  // center of mass and inertia tensor for the entire mesh
  void getInertiaParameters(double & mass, Vec3d & centerOfMass, Mat3d & inertiaTensor);

  // center is the geometric center of all vertices; radius is tightest fitting
  // sphere centered at the center
  void getMeshGeometricParameters(double * centerX, double * centerY, double * centerZ, double * radius);

  // proximity queries

  // finds the closest element to the given position (using linear scan);
  // distance to a element is defined as distance to its center
  int getClosestElement(Vec3d pos); 

  // finds the closest vertex to the given position (using linear scan)
  int getClosestVertex(Vec3d pos); 

  // finds the element that containts the given position (using linear scan); if
  // such element does not exist, -1 is returned
  int getContainingElement(Vec3d pos); 

  // true if given element contain given position, false otherwise
  virtual bool containsVertex(int element, Vec3d pos) = 0; 

  // exports the mesh geometric information (for external usage)
  void dumpMeshData(int * numVertices, double ** vertices, int * numElements, int * numElementVertices, int ** elements);

  // computes the gravity vector (different forces on different mesh vertices
  // due to potentially varying mass densities) gravityForce must be a
  // pre-allocated vector of length 3xnumVertices()
  void computeGravity(double * gravityForce, double g=9.81, bool addForce=false);

  // edge queries
  virtual int numElementEdges()const = 0;

  // edgeBuffer must be pre-allocated, of size 2 x numElementEdges()
  virtual void getElementEdges(int el, int * edgeBuffer)const = 0;

  // computes the mass matrix of a single element
  // note: to compute the mass matrix for the entire mesh, use generateMassMatrix.h
  // massMatrix is numElementVertices_ x numElementVertices_
  virtual void computeElementMassMatrix(int element, double * massMatrix) = 0; 

  // (permanently) applies the deformation to the vertices of the mesh
  void applyDeformation(double * u);

  // === interpolation ===
  // the interpolant is a triple (numTargetLocations, vertices, weights)

  // Builds the interpolant
  // input is a list of 3D target locations where the interpolant will be built
  // e.g. those could be vertices of a triangle mesh embedded into the sqmesh
  // each location is a 3-vector, i.e., 3 consecutive double-precision vectors
  // if zeroThreshold is set positive, than for any target location that is 
  //   more than zeroThreshold away from the closest voxel, 
  //   all weights will be set to zero; this is useful, e.g. to 
  //   stabilize locations far away from your mesh
  // output: vertices and weights arrays
  // vertices: gives a list of integer indicies of the vertices of the element
  //   closest to the target location (numElementVertices entries per target location, one for each element vertex)
  //   note: if target location is inside a voxel, that voxel will be chosen as closest
  // weights: a list of numElementVertices_ weights, as per the numElementVertices_ vertices of each element (weights sum to 1)

  // returns the number of targetLocations that were outside of the mesh
  int buildInterpolationWeights(int numTargetLocations, const float * targetLocations,
								int ** vertices, double ** weights, double zeroThreshold = -1.0, 
								std::vector<int> * elementList = NULL, int verbose=0);

  // looks at the first line of "filename" to determine "numElementVertices" for
  // this particular interpolant
  static int getNumInterpolationElementVertices(char * filename); 

  // returns 0 on success
  static int loadInterpolationWeights(const char * filename, int numTargetLocations,int numElementVertices, 
									  int ** vertices, double ** weights);
  static int saveInterpolationWeights(const char * filename, int numTargetLocations,int numElementVertices, 
									  const int * vertices, const double * weights);

  // interpolates 3D vector data from vertices of the volumetric mesh (data
  // given in u) to the target locations (output goes into uTarget) e.g., use
  // this to interpolate deformation from the volumetric mesh to a triangle mesh
  static void interpolate(const double * u, float * uTarget, 
						  int numTargetLocations, int numElementVertices, 
						  const int * vertices, const double * weights);

  // computes barycentric weights of the given position with respect to the given element
  virtual void computeBarycentricWeights(int element, Vec3d pos, double * weights) = 0;

protected:
  int numElementVertices_;

  int numVertices_;
  Vec3d ** vertices_;

  int numElements_;
  MeshElement ** elements_;

  int numMaterials_;
  int numElementSets_;
  int numSolidSections_;
  Material ** materials_;
  ElementSet ** elementSets_; // elements inside these sets are 1-indexed
  SolidSection ** solidSections_;

  int save(char * filename, char * elementName);

  // parses in the mesh, and returns the string corresponding to the element type
  // (pointer is allocated inside the routine, or set to NULL upon failure)
  VolumetricMesh(const char * filename, int numElementVertices, char ** elementType);
  VolumetricMesh(int numElementVertices) { numElementVertices_ = numElementVertices;}
  void PropagateMaterialsToElements();

  char * temp; // auxiliary
};

typedef boost::shared_ptr< VolumetricMesh > pVolumetricMesh; 
typedef boost::shared_ptr< const VolumetricMesh > pVolumetricMesh_const;

#endif

