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

#include <float.h>
#include <string.h>
#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>
#include "sqmesh.h"
#include "matrixMacros.h"
using namespace std;

const char SQMesh::elementType_[96]="C3D8";

SQMesh::SQMesh(const char * filename) : VolumetricMesh(filename, 8, &temp)
{
  if (strcmp(temp, elementType_) != 0)
  {
    printf("Error: unknown element type %s encountered.\n", temp);
    throw 11;
  } 
  free(temp);

  // set cube size
  cubeSize_ = len(*vertex(0,1)-*vertex(0,0));
}

SQMesh::SQMesh(int numVertices, double * vertices,
               int numElements, int numElementVertices, int * elements,
               double E, double nu, double density): VolumetricMesh(numVertices, vertices, numElements, numElementVertices, elements, E, nu, density)
{
  if (numElementVertices_ != 8)
    throw 1;

  cubeSize_ = len(*vertex(0,1)-*vertex(0,0));
}

SQMesh::SQMesh(SQMesh & sqmesh, int numElements, int * elements, map<int,int> * vertexMap_): VolumetricMesh(sqmesh, numElements, elements, vertexMap_)
{
  cubeSize_ = sqmesh.cubeSize();
}

SQMesh::~SQMesh() {}

int SQMesh::save(char * filename)
{
  return VolumetricMesh::save(filename, (char*)elementType_);
}

bool SQMesh::containsVertex(int element, Vec3d pos)
{
  Vec3d bmin = *(vertex(element, 0));
  Vec3d bmax = *(vertex(element, 6));

  return (pos[0] >= bmin[0]) && (pos[1] >= bmin[1]) && (pos[2] >= bmin[2]) &&
         (pos[0] <= bmax[0]) && (pos[1] <= bmax[1]) && (pos[2] <= bmax[2]);
}

// interpolates given sq mesh vertex 3D data to the destination locations
// vertexData size must a vector of length 3 * nElements
// destMatrix must be a vector of length 3 * numLocations
// if a given location does not belong to any element, extrapolation to the nearest element will be used
// vertices more than distanceThreshold away from any element vertex are assigned zero data
int SQMesh::interpolateData(double * vertexData, int numInterpolationLocations, int r, double * interpolationLocations, double * destMatrix, double distanceThreshold)
{
  assert(numElementVertices_ == 8);

  int numExternalVertices = 0;

  for (int i=0; i < numInterpolationLocations; i++) // over all interpolation locations
  {
    if (i % 100 == 0)
    {
      printf("%d ", i); fflush(NULL);
    }

    Vec3d pos = Vec3d(interpolationLocations[3*i+0], interpolationLocations[3*i+1], interpolationLocations[3*i+2]);

    // find element containing pos
    int element = getClosestElement(pos);

    if (!containsVertex(element, pos))
      numExternalVertices++;

    if (distanceThreshold > 0)
    {
      // check whether vertex is close enough to the sqmesh
      double minDistance = DBL_MAX;
      int numElementVertices_ = numElementVertices();
      int assignedZero = 0;
      for(int ii=0; ii< numElementVertices_; ii++)
      {
        Vec3d * vpos = vertex(element, ii);
        if (len(*vpos-pos) < minDistance)
        {
          minDistance = len(*vpos-pos);
        }
      }

      if (minDistance > distanceThreshold)
      {
        // assign zero data
        for(int j=0; j<r; j++)
        {
          destMatrix[ELT(3*numInterpolationLocations,3*i+0,j)] = 0;
          destMatrix[ELT(3*numInterpolationLocations,3*i+1,j)] = 0;
          destMatrix[ELT(3*numInterpolationLocations,3*i+2,j)] = 0;
        }
        assignedZero++;
        continue;
      }
    }

    // compute barycentric coordinates
    Vec3d w = pos - *(vertex(element, 0));
    double alpha = w[0] / cubeSize_;
    double beta = w[1] / cubeSize_;
    double gamma = w[2] / cubeSize_;

    double f000 = (1-alpha)*(1-beta)*(1-gamma);
    double f100 = (alpha)*(1-beta)*(1-gamma);
    double f110 = (alpha)*(beta)*(1-gamma);
    double f010 = (1-alpha)*(beta)*(1-gamma);

    double f001 = (1-alpha)*(1-beta)*(gamma);
    double f101 = (alpha)*(1-beta)*(gamma);
    double f111 = (alpha)*(beta)*(gamma);
    double f011 = (1-alpha)*(beta)*(gamma);

    int v000 = vertexIndex(element, 0);
    int v100 = vertexIndex(element, 1);
    int v110 = vertexIndex(element, 2);
    int v010 = vertexIndex(element, 3);
    int v001 = vertexIndex(element, 4);
    int v101 = vertexIndex(element, 5);
    int v111 = vertexIndex(element, 6);
    int v011 = vertexIndex(element, 7);

    for(int j=0; j<r; j++)
    {
      Vec3d data000 = Vec3d (vertexData[ELT(3*numVertices_, 3 * v000 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v000 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v000 + 2, j)]);
      Vec3d data100 = Vec3d (vertexData[ELT(3*numVertices_, 3 * v100 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v100 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v100 + 2, j)]);
      Vec3d data110 = Vec3d (vertexData[ELT(3*numVertices_, 3 * v110 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v110 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v110 + 2, j)]);
      Vec3d data010 = Vec3d (vertexData[ELT(3*numVertices_, 3 * v010 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v010 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v010 + 2, j)]);

      Vec3d data001 = Vec3d (vertexData[ELT(3*numVertices_, 3 * v001 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v001 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v001 + 2, j)]);
      Vec3d data101 = Vec3d (vertexData[ELT(3*numVertices_, 3 * v101 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v101 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v101 + 2, j)]);
      Vec3d data111 = Vec3d (vertexData[ELT(3*numVertices_, 3 * v111 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v111 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v111 + 2, j)]);
      Vec3d data011 = Vec3d (vertexData[ELT(3*numVertices_, 3 * v011 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v011 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v011 + 2, j)]);

      Vec3d interpolatedData = f000 * data000 + f100 * data100 + f110 * data110 + f010 * data010 +
                               f001 * data001 + f101 * data101 + f111 * data111 + f011 * data011 ;

      destMatrix[ELT(3*numInterpolationLocations,3*i+0,j)] = interpolatedData[0];
      destMatrix[ELT(3*numInterpolationLocations,3*i+1,j)] = interpolatedData[1];
      destMatrix[ELT(3*numInterpolationLocations,3*i+2,j)] = interpolatedData[2];
    }
  }

  printf("\n");

  return numExternalVertices;
}

// computes approximation to the normal correction for the given deformation
// vertexData size must a vector of length 3 * nElements
// normalCorrection (output) must be a vector of length 3 * numLocations
// staticNormals must be a vector of length 3 * numLocations
// vertices more than distanceThreshold away from any element vertex are assigned zero data
int SQMesh::normalCorrection(double * vertexData, int numInterpolationLocations, int r, double * interpolationLocations, double * staticNormals, double * normalCorrection, double distanceThreshold)
{
  assert(numElementVertices_ == 8);

  int numExternalVertices = 0;

  for (int i=0; i < numInterpolationLocations; i++) // over all interpolation locations
  {

    if (i % 100 == 0)
    {
      printf("%d ", i); fflush(NULL);
    }

    Vec3d pos = Vec3d(interpolationLocations[3*i+0], interpolationLocations[3*i+1], interpolationLocations[3*i+2]);

    // find element containing pos
    int element = getClosestElement(pos);

    if (!containsVertex(element, pos))
      numExternalVertices++;

    if (distanceThreshold > 0)
    {
      // check whether vertex is close enough to the sqmesh
      double minDistance = DBL_MAX;
      int numElementVertices_ = numElementVertices();
      int assignedZero = 0;
      for(int ii=0; ii< numElementVertices_; ii++)
      {
        Vec3d * vpos = vertex(element, ii);
        if (len(*vpos-pos) < minDistance)
        {
          minDistance = len(*vpos-pos);
        }
      }

      if (minDistance > distanceThreshold)
      {
        // assign zero data
        for(int j=0; j<r; j++)
        {
          normalCorrection[ELT(3*numInterpolationLocations,3*i+0,j)] = 0;
          normalCorrection[ELT(3*numInterpolationLocations,3*i+1,j)] = 0;
          normalCorrection[ELT(3*numInterpolationLocations,3*i+2,j)] = 0;
        }
        assignedZero++;

        continue;
      }
    }

    // compute barycentric coordinates
    Vec3d w = pos - *(vertex(element, 0));
    double alpha = w[0] / cubeSize_;
    double beta = w[1] / cubeSize_;
    double gamma = w[2] / cubeSize_;

    //double f000 = (1-alpha)*(1-beta)*(1-gamma);
    //double f100 = (alpha)*(1-beta)*(1-gamma);
    //double f110 = (alpha)*(beta)*(1-gamma);
    //double f010 = (1-alpha)*(beta)*(1-gamma);

    //double f001 = (1-alpha)*(1-beta)*(gamma);
    //double f101 = (alpha)*(1-beta)*(gamma);
    //double f111 = (alpha)*(beta)*(gamma);
    //double f011 = (1-alpha)*(beta)*(gamma);

    int v000 = vertexIndex(element, 0);
    int v100 = vertexIndex(element, 1);
    int v110 = vertexIndex(element, 2);
    int v010 = vertexIndex(element, 3);
    int v001 = vertexIndex(element, 4);
    int v101 = vertexIndex(element, 5);
    int v111 = vertexIndex(element, 6);
    int v011 = vertexIndex(element, 7);

    Vec3d gradf000(1.0 / cubeSize_ * -(1-beta)*(1-gamma), 1.0 / cubeSize_ * -(1-alpha)*(1-gamma), 1.0 / cubeSize_ * -(1-alpha)*(1-beta));
    Vec3d gradf100(1.0 / cubeSize_ * (1-beta)*(1-gamma), 1.0 / cubeSize_ * -alpha*(1-gamma), 1.0 / cubeSize_ * -alpha*(1-beta));
    Vec3d gradf110(1.0 / cubeSize_ * beta*(1-gamma), 1.0 / cubeSize_ * alpha*(1-gamma), 1.0 / cubeSize_ * -alpha*beta);
    Vec3d gradf010(1.0 / cubeSize_ * -beta*(1-gamma), 1.0 / cubeSize_ * (1-alpha)*(1-gamma), 1.0 / cubeSize_ * (1-alpha)*-beta);

    Vec3d gradf001(1.0 / cubeSize_ * -(1-beta)*gamma, 1.0 / cubeSize_ * -(1-alpha)*gamma, 1.0 / cubeSize_ * (1-alpha)*(1-beta));
    Vec3d gradf101(1.0 / cubeSize_ * (1-beta)*gamma, 1.0 / cubeSize_ * -alpha*gamma, 1.0 / cubeSize_ * alpha*(1-beta));
    Vec3d gradf111(1.0 / cubeSize_ * beta*gamma, 1.0 / cubeSize_ * alpha*gamma, 1.0 / cubeSize_ * alpha*beta);
    Vec3d gradf011(1.0 / cubeSize_ * -beta*gamma, 1.0 / cubeSize_ * (1-alpha)*gamma, 1.0 / cubeSize_ * (1-alpha)*beta);

    Vec3d normal = Vec3d(staticNormals[3*i+0], staticNormals[3*i+1], staticNormals[3*i+2]);

    for(int j=0; j<r; j++)
    {
      Vec3d u000 = Vec3d(vertexData[ELT(3*numVertices_, 3 * v000 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v000 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v000 + 2, j)]);
      Vec3d u100 = Vec3d(vertexData[ELT(3*numVertices_, 3 * v100 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v100 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v100 + 2, j)]);
      Vec3d u110 = Vec3d(vertexData[ELT(3*numVertices_, 3 * v110 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v110 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v110 + 2, j)]);
      Vec3d u010 = Vec3d(vertexData[ELT(3*numVertices_, 3 * v010 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v010 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v010 + 2, j)]);

      Vec3d u001 = Vec3d(vertexData[ELT(3*numVertices_, 3 * v001 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v001 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v001 + 2, j)]);
      Vec3d u101 = Vec3d(vertexData[ELT(3*numVertices_, 3 * v101 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v101 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v101 + 2, j)]);
      Vec3d u111 = Vec3d(vertexData[ELT(3*numVertices_, 3 * v111 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v111 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v111 + 2, j)]);
      Vec3d u011 = Vec3d(vertexData[ELT(3*numVertices_, 3 * v011 + 0, j)], vertexData[ELT(3*numVertices_, 3 * v011 + 1, j)], vertexData[ELT(3*numVertices_, 3 * v011 + 2, j)]);

      Vec3d coef(0,0,0);
      coef += dot(gradf000,normal) * u000;
      coef += dot(gradf100,normal) * u100;
      coef += dot(gradf110,normal) * u110;
      coef += dot(gradf010,normal) * u010;
      coef += dot(gradf001,normal) * u001;
      coef += dot(gradf101,normal) * u101;
      coef += dot(gradf111,normal) * u111;
      coef += dot(gradf011,normal) * u011;

      normalCorrection[ELT(3*numInterpolationLocations,3*i+0,j)] = coef[0];
      normalCorrection[ELT(3*numInterpolationLocations,3*i+1,j)] = coef[1];
      normalCorrection[ELT(3*numInterpolationLocations,3*i+2,j)] = coef[2];
    }
  }

  printf("\n");

  return numExternalVertices;
}

int SQMesh::deformationGradientModes(double * sqMeshVertexModes, int r, Vec3d pos, double * F)
{
  // find element containing pos
  int element = getClosestElement(pos);

  int externalVertex = !containsVertex(element, pos);

  // compute barycentric coordinates
  Vec3d w = pos - *(vertex(element, 0));
  double alpha = w[0] / cubeSize_;
  double beta = w[1] / cubeSize_;
  double gamma = w[2] / cubeSize_;

  //double f000 = (1-alpha)*(1-beta)*(1-gamma);
  //double f100 = (alpha)*(1-beta)*(1-gamma);
  //double f110 = (alpha)*(beta)*(1-gamma);
  //double f010 = (1-alpha)*(beta)*(1-gamma);

  //double f001 = (1-alpha)*(1-beta)*(gamma);
  //double f101 = (alpha)*(1-beta)*(gamma);
  //double f111 = (alpha)*(beta)*(gamma);
  //double f011 = (1-alpha)*(beta)*(gamma);

  int v000 = vertexIndex(element, 0);
  int v100 = vertexIndex(element, 1);
  int v110 = vertexIndex(element, 2);
  int v010 = vertexIndex(element, 3);
  int v001 = vertexIndex(element, 4);
  int v101 = vertexIndex(element, 5);
  int v111 = vertexIndex(element, 6);
  int v011 = vertexIndex(element, 7);

  Vec3d gradf000(1.0 / cubeSize_ * -(1-beta)*(1-gamma), 1.0 / cubeSize_ * -(1-alpha)*(1-gamma), 1.0 / cubeSize_ * -(1-alpha)*(1-beta));
  Vec3d gradf100(1.0 / cubeSize_ * (1-beta)*(1-gamma), 1.0 / cubeSize_ * -alpha*(1-gamma), 1.0 / cubeSize_ * -alpha*(1-beta));
  Vec3d gradf110(1.0 / cubeSize_ * beta*(1-gamma), 1.0 / cubeSize_ * alpha*(1-gamma), 1.0 / cubeSize_ * -alpha*beta);
  Vec3d gradf010(1.0 / cubeSize_ * -beta*(1-gamma), 1.0 / cubeSize_ * (1-alpha)*(1-gamma), 1.0 / cubeSize_ * (1-alpha)*-beta);

  Vec3d gradf001(1.0 / cubeSize_ * -(1-beta)*gamma, 1.0 / cubeSize_ * -(1-alpha)*gamma, 1.0 / cubeSize_ * (1-alpha)*(1-beta));
  Vec3d gradf101(1.0 / cubeSize_ * (1-beta)*gamma, 1.0 / cubeSize_ * -alpha*gamma, 1.0 / cubeSize_ * alpha*(1-beta));
  Vec3d gradf111(1.0 / cubeSize_ * beta*gamma, 1.0 / cubeSize_ * alpha*gamma, 1.0 / cubeSize_ * alpha*beta);
  Vec3d gradf011(1.0 / cubeSize_ * -beta*gamma, 1.0 / cubeSize_ * (1-alpha)*gamma, 1.0 / cubeSize_ * (1-alpha)*beta);

  for(int j=0; j<r; j++)
  {
    Vec3d u000 = Vec3d(sqMeshVertexModes[ELT(3*numVertices_, 3 * v000 + 0, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v000 + 1, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v000 + 2, j)]);
    Vec3d u100 = Vec3d(sqMeshVertexModes[ELT(3*numVertices_, 3 * v100 + 0, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v100 + 1, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v100 + 2, j)]);
    Vec3d u110 = Vec3d(sqMeshVertexModes[ELT(3*numVertices_, 3 * v110 + 0, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v110 + 1, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v110 + 2, j)]);
    Vec3d u010 = Vec3d(sqMeshVertexModes[ELT(3*numVertices_, 3 * v010 + 0, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v010 + 1, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v010 + 2, j)]);

    Vec3d u001 = Vec3d(sqMeshVertexModes[ELT(3*numVertices_, 3 * v001 + 0, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v001 + 1, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v001 + 2, j)]);
    Vec3d u101 = Vec3d(sqMeshVertexModes[ELT(3*numVertices_, 3 * v101 + 0, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v101 + 1, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v101 + 2, j)]);
    Vec3d u111 = Vec3d(sqMeshVertexModes[ELT(3*numVertices_, 3 * v111 + 0, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v111 + 1, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v111 + 2, j)]);
    Vec3d u011 = Vec3d(sqMeshVertexModes[ELT(3*numVertices_, 3 * v011 + 0, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v011 + 1, j)], sqMeshVertexModes[ELT(3*numVertices_, 3 * v011 + 2, j)]);

    Mat3d FMode(0.0);
    FMode += tensorProduct(u000, gradf000);
    FMode += tensorProduct(u100, gradf100);
    FMode += tensorProduct(u110, gradf110);
    FMode += tensorProduct(u010, gradf010);
    FMode += tensorProduct(u001, gradf001);
    FMode += tensorProduct(u101, gradf101);
    FMode += tensorProduct(u111, gradf111);
    FMode += tensorProduct(u011, gradf011);

    for(int ii=0; ii<3; ii++)
      for(int jj=0; jj<3; jj++)
      {
        // 9 x r matrix
        F[ELT(9, 3*ii+jj, j)] = FMode[ii][jj];
      }
  }

  return externalVertex;
}

void SQMesh::computeBarycentricWeights(int el, Vec3d pos, double * weights)
{
  // compute barycentric coordinates
  Vec3d w = pos - *(vertex(el, 0));
  double alpha = w[0] / cubeSize_;
  double beta = w[1] / cubeSize_;
  double gamma = w[2] / cubeSize_;

  weights[0] = (1-alpha)*(1-beta)*(1-gamma); // f000
  weights[1] = (alpha)*(1-beta)*(1-gamma); // f100
  weights[2] = (alpha)*(beta)*(1-gamma); // f110
  weights[3] = (1-alpha)*(beta)*(1-gamma); // f010

  weights[4] = (1-alpha)*(1-beta)*(gamma); // f001
  weights[5] = (alpha)*(1-beta)*(gamma); // f101
  weights[6] = (alpha)*(beta)*(gamma); // f111
  weights[7] = (1-alpha)*(beta)*(gamma); // f011
}

double SQMesh::getElementVolume(int el)const{

  return cubeSize_ * cubeSize_ * cubeSize_;
}

void SQMesh::getElementInertiaTensor(int el, Mat3d & inertiaTensor)
{
  // inertia tensor of a cube with side a, and mass m is spherical, with
  // diagonal elements 1.0 / 6 * m * a * a
  double cubeSize5 = cubeSize_ * cubeSize_ * cubeSize_ * cubeSize_ * cubeSize_;

  inertiaTensor = cubeSize5 * Mat3d( 1.0/6 , 0 , 0, 0 , 1.0/6 , 0, 0 , 0 , 1.0/6 );
}

int SQMesh::numElementEdges()const
{
  return 12;
}

void SQMesh::getElementEdges(int el, int * edgeBuffer)const
{
  int v[8];
  for(int i=0; i<8; i++)
    v[i] = vertexIndex(el,i);

  int edgeMask[12][2] = {
   { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 },
   { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 },
   { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } } ;

  for(int edge=0; edge<12; edge++)
  {
    edgeBuffer[2*edge+0] = v[edgeMask[edge][0]];
    edgeBuffer[2*edge+1] = v[edgeMask[edge][1]];
  }
}

void SQMesh::computeElementMassMatrix(int el, double * massMatrix)
{
  double singleMassMatrix[64] =
  { 0.0370370335876942, 0.0185185167938471, 0.00925925839692354, 0.0185185167938471,
    0.0185185167938471, 0.00925925839692354, 0.00462962919846177, 0.00925925839692354,
    0.0185185167938471, 0.0370370335876942, 0.0185185167938471, 0.00925925839692354,
    0.00925925839692354, 0.0185185167938471, 0.00925925839692354, 0.00462962919846177,
    0.00925925839692354, 0.0185185167938471, 0.0370370335876942, 0.0185185167938471,
    0.00462962919846177, 0.00925925839692354, 0.0185185167938471, 0.00925925839692354,
    0.0185185167938471, 0.00925925839692354, 0.0185185167938471, 0.0370370335876942,
    0.00925925839692354, 0.00462962919846177, 0.00925925839692354, 0.0185185167938471,
    0.0185185167938471, 0.00925925839692354, 0.00462962919846177, 0.00925925839692354,
    0.0370370335876942, 0.0185185167938471, 0.00925925839692354, 0.0185185167938471,
    0.00925925839692354, 0.0185185167938471, 0.00925925839692354, 0.00462962919846177,
    0.0185185167938471, 0.0370370335876942, 0.0185185167938471, 0.00925925839692354,
    0.00462962919846177, 0.00925925839692354, 0.0185185167938471, 0.00925925839692354,
    0.00925925839692354, 0.0185185167938471, 0.0370370335876942, 0.0185185167938471,
    0.00925925839692354, 0.00462962919846177, 0.00925925839692354, 0.0185185167938471,
    0.0185185167938471, 0.00925925839692354, 0.0185185167938471, 0.0370370335876942
  };

  double density = element(el)->density();
  double voxelSpacing = cubeSize();
  double voxelSpacing3 = voxelSpacing * voxelSpacing * voxelSpacing;
  double factor = density * voxelSpacing3;

  for(int i=0; i<64; i++)
    massMatrix[i] = factor * singleMassMatrix[i];
}

