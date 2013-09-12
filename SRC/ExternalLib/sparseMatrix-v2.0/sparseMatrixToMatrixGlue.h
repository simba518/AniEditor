/*************************************************************************
 *                                                                       *
 * "sparseMatrix" library , Copyright (C) 2007 CMU, 2009 MIT             *
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

#ifndef _SPARSEMATRIXTOMATRIXGLUE_H_
#define _SPARSEMATRIXTOMATRIXGLUE_H_

#include "sparseMatrix.h"

/*
  Export a sparse matrix into a dense matrix of the "Matrix" class ("matrix" library).
  Note: this routine requires the "matrix" library. If you only need sparse matrices,
  you don't need to use this file.
*/

const BARBIC_LIB::Matrix<double>operator*(const BARBIC_LIB::SparseMatrix & A , 
										  const BARBIC_LIB::Matrix<double> & mtx);

#endif
