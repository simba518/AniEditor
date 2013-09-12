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

#include "volumetricMeshLoader.h"
#include "sqmesh.h"
#include "tetmesh.h"

VolumetricMesh * VolumetricMeshLoader::load(const char * filename)
{
  char * elementType;
  VolumetricMesh * volumetricMesh = NULL;

  if (VolumetricMesh::getElementType(filename, &elementType) == NULL)
    return NULL;

  if (strcmp(elementType,SQMesh::elementType()) == 0)
    volumetricMesh = new SQMesh(filename); 

  if (strcmp(elementType,TetMesh::elementType()) == 0)
    volumetricMesh = new TetMesh(filename); 

  free(elementType);

  return volumetricMesh;
}

