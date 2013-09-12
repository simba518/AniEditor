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

/*

  This class is a container for a cube volumetric 3D mesh. 
  See also volumetricMesh.h. 

  All cubes are of equal size, aligned with coordinate system axes, and 
  follow a grid pattern. They are typically obtained by voxelizing a 
  given input triangle geometry.

  To generate a sqmesh from an input triangle mesh (optionally flood-filling 
  interior chambers), you can use the "Large Modal Deformation Factory" 
  application ( http://www.jernejbarbic.com/code ) .

*/

#ifndef _SQMESH_H_
#define _SQMESH_H_

#include "volumetricMesh.h"

// see also volumetricMesh.h for a description of the routines

class SQMesh : public VolumetricMesh
{
public:
  
  // loads the mesh from a text file (note: input format is similar to that of Abaqus)
  SQMesh(const char * filename); 

  // constructs a mesh from the given vertices and elements, with a single solid section and material
  SQMesh(int numVertices, double * vertices,
         int numElements, int numElementVertices, int * elements,
         double E, double nu, double density);

  // creates a mesh consisting of the specified element subset of the given sqmesh
  SQMesh(SQMesh & mesh, int numElements, int * elements, std::map<int,int> * vertexMap = NULL);

  virtual ~SQMesh();

  // saves the mesh to a text file (a format similar to Abaqus is used; only geometry is saved; not material information)
  virtual int save(char * filename);

  // misc queries
  static const char * elementType() { return elementType_; }
  virtual const char * getElementType() { return elementType(); }
  inline double cubeSize() { return cubeSize_; }
  virtual double getElementVolume(int el)const;
  virtual void getElementInertiaTensor(int el, Mat3d & inertiaTensor);
  virtual bool containsVertex(int element, Vec3d pos); // true if given element contain given position, false otherwise

  // edge queries
  virtual int numElementEdges()const;  
  virtual void getElementEdges(int el, int * edgeBuffer)const;

  // === interpolation ===
  virtual void computeBarycentricWeights(int el, Vec3d pos, double * weights);
  virtual void computeElementMassMatrix(int el, double * massMatrix);

  int interpolateData(double * volumetricMeshVertexData, int numLocations, int r, double * interpolationLocations, double * destMatrix, double zeroThreshold = -1.0);

  // computes approximation to the normal correction for the given deformations
  // vertexData size must a matrix of size 3 * nElements x r
  // destMatrix must be a vector of length 3 * numLocations x r
  // staticNormals must be a vector of length 3 * numLocations
  // returns the number of vertices that were not contained inside any element
  // vertices more than distanceThreshold away from any element vertex are assigned zero data
  int normalCorrection(double * vertexData, int numLocations, int r, double * interpolationLocations, double * staticNormals, double * normalCorrection, double zeroThreshold = -1.0); // note: this routine could be promoted to volumetricMesh.h
  int deformationGradientModes(double * vertexData, int r, Vec3d pos, double * F);

protected:
  double cubeSize_;

  static const char elementType_[96];
};

#endif

