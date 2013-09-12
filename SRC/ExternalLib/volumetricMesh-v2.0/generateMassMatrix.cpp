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

#include "generateMassMatrix.h"
using namespace BARBIC_LIB;

void GenerateMassMatrix::computeMassMatrix(
  VolumetricMesh * volumetricMesh, SparseMatrix ** massMatrix, bool inflate3Dim)
{
  int n = volumetricMesh->numVertices();
  int numElementVertices = volumetricMesh->numElementVertices();
  double * buffer = (double*) malloc (sizeof(double) * numElementVertices * numElementVertices);

  SparseMatrixOutline * massMatrixOutline;
  if (!inflate3Dim)
  {
    massMatrixOutline = new SparseMatrixOutline(n);
    for(int el=0; el <volumetricMesh->numElements(); el++)
    {
      volumetricMesh->computeElementMassMatrix(el, buffer);
      for(int i=0; i < numElementVertices; i++)
        for(int j=0; j < numElementVertices; j++)
        {
          massMatrixOutline->AddEntry(volumetricMesh->vertexIndex(el,i),volumetricMesh->vertexIndex(el,j),  buffer[numElementVertices * j + i]);
        }
    }
  }
  else
  {
    massMatrixOutline = new SparseMatrixOutline(3*n);
    for(int el=0; el <volumetricMesh->numElements(); el++)
    {
      volumetricMesh->computeElementMassMatrix(el, buffer);
      for(int i=0; i < numElementVertices; i++)
        for(int j=0; j < numElementVertices; j++)
        {
          double entry = buffer[numElementVertices * j + i];
          int indexi = volumetricMesh->vertexIndex(el,i);
          int indexj = volumetricMesh->vertexIndex(el,j);
          massMatrixOutline->AddEntry(3*indexi+0, 3*indexj+0, entry);
          massMatrixOutline->AddEntry(3*indexi+1, 3*indexj+1, entry);
          massMatrixOutline->AddEntry(3*indexi+2, 3*indexj+2, entry);
        }
    }
  }

  (*massMatrix) = new SparseMatrix(massMatrixOutline);
  delete(massMatrixOutline);

  free(buffer);
}

