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

#include "meshElement.h"
#include <stdlib.h>

MeshElement::MeshElement(double E, double nu, double density, int numVertices, int * vertices):
      E_(E), nu_(nu), density_(density), numVertices_(numVertices)
{ 
  vertices_ = (int*) malloc (sizeof(int) * numVertices_);
  for (int i=0; i<numVertices_; i++)
    vertices_[i] = vertices[i];
}

MeshElement::MeshElement(MeshElement & meshElement)
{
  E_ = meshElement.E();
  nu_ = meshElement.nu();
  density_ = meshElement.density();
  numVertices_ = meshElement.numVertices_;
  vertices_ = (int*) realloc(vertices_, sizeof(int) * numVertices_);
  for(int i=0; i<numVertices_; i++)
    vertices_[i] = meshElement.vertex(i);
}

MeshElement::~MeshElement()
{
  free(vertices_);
}

