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
  This class is a container for a tetrahedral volumetric 3D mesh. See 
  also volumetricMesh.h. The tetrahedra can take arbitrary shapes (not 
  limited to only a few shapes).
*/

#ifndef _TETMESH_H_
#define _TETMESH_H_

#include "volumetricMesh.h"

// see also volumetricMesh.h for a description of the routines

class TetMesh : public VolumetricMesh
{
public:
  // loads the mesh from a text file 
  // (input formut is similar to Abaqus, see the provided example)
  TetMesh(const char * filename);

  // loads a file of a "special" (non-Abaqus-related) type
  // currently one such special format is supported:
  // specialFileType=0: 
  //   the ".ele" and ".node" format, used by TetGen, 
  //   "filename" is the basename, e.g., passing "mesh" will load the mesh from "mesh.ele" and "mesh.node" (default material parameters will be used)
  TetMesh(char * filename, int specialFileType); 

  virtual ~TetMesh();

  virtual int save(char * filename);

  static const char * elementType() { return elementType_; }
  virtual const char * getElementType() { return elementType(); }
  virtual double getElementVolume(int el)const;
  virtual void getElementInertiaTensor(int el, Mat3d & inertiaTensor);
  virtual bool containsVertex(int element, Vec3d pos); // true if given element contain given position, false otherwise
  virtual void computeBarycentricWeights(int el, Vec3d pos, double * weights);
  virtual void computeElementMassMatrix(int element, double * massMatrix);

  static double determinantHelper(Vec3d a, Vec3d b, Vec3d c, Vec3d d);

  virtual int numElementEdges()const;
  virtual void getElementEdges(int el, int * edgeBuffer)const;

protected:

  void computeElementMassMatrixHelper(Vec3d a, Vec3d b, Vec3d c, Vec3d d, double * buffer);

  static const char elementType_[96];
};

#endif

