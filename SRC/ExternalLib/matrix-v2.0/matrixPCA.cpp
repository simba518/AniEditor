/*************************************************************************
 *                                                                       *
 * "matrix" library , Copyright (C) 2007 CMU, 2009 MIT                   *
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
#include <log4cplus/logger.h>
#include <log4cplus/loggingmacros.h>
#include <loggingmacrosext.h>
using namespace log4cplus;

#include "matrixPCA.h"
#include "matrixMacros.h"
#include "matrixBLAS.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <lapack-headers.h>

void DoTresholding_Epsilon(double * singularValues, int numberOfSingularValues, int * r, double epsilon)
{
  // r = number of retained components
  // determine r

  // first-singular value-based tresholding
  double energy = 0;
  *r = 0;
  while ((singularValues[*r] >= epsilon * singularValues[0]) && (*r < numberOfSingularValues) )
  {
    energy += singularValues[*r] * singularValues[*r];
    (*r)++;
  }

  // determine total energy
  double totalEnergy = energy;
  int s = *r;
  while (s < numberOfSingularValues)
  {
    totalEnergy += singularValues[s] * singularValues[s];
    s++;
  }

  printf ("Retained components: %d , PCA epsilon treshold = %f , relative energy retained = %f .\n",*r,epsilon,energy/totalEnergy);
}

void DoTresholding_NumberOfModes(double * singularValues, int numberOfSingularValues, int * r, int rDesired)
{
  // r = number of retained components
  // determine r

  // first-singular value-based tresholding
  double energy = 0;
  *r = 0;
  while ((*r < rDesired) && (*r < numberOfSingularValues) )
  {
    energy += singularValues[*r] * singularValues[*r];
    (*r)++;
  }

  // determine total energy
  double totalEnergy = energy;
  int s = *r;
  while (s < numberOfSingularValues)
  {
    totalEnergy += singularValues[s] * singularValues[s];
    s++;
  }
  Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("ExternalLib"));
  LOG4CPLUS_DEBUG(logger,"Retained components = " << *r);
  LOG4CPLUS_DEBUG(logger,"requested components = " << rDesired);
  LOG4CPLUS_DEBUG(logger,"relative energy retained = " << energy/totalEnergy);
}

int MatrixPCA(ThresholdingSpecification * thresholdingSpecification,
              int m, int n, double * A, int * r){

  Logger logger = Logger::getInstance(LOG4CPLUS_TEXT("ExternalLib"));
  if (!A)
  {
	LOG4CPLUS_ERROR(logger,"input matrix is NULL.");
    return -1;
  }

  bool transpose = false;
  if (m > n)
  {
    transpose = true;
    InPlaceTransposeMatrix(m,n,A);

    // swap m,n
    int bufferi = m;
    m = n;
    n = bufferi;
  }

  // do SVD

  char jobu  = 'O';//overwrites A with U (left singular vectors)
  //char jobu  = 'N';
  char jobvt = 'S';//all rows returned in VT
  //char jobvt = 'N';

  int ldA = m;
  int ldU = m;
  int lwork = 64*MAX(3*MIN( m, n)+MAX(m,n), 5*MIN(m,n)-4);
  double * work = (double*) malloc (sizeof(double) * lwork);
  if (!work)
  {
    LOG4CPLUS_ERROR(logger,"failed to allocate workspace.");
    return -2;
  }
  LOG4CPLUS_DEBUG(logger, "Workspace size is: " << 1.0 * lwork * sizeof(int) / 1024 / 1024 <<"Mb");

  // allocate array for singular vectors
  double * S = (double *) malloc (sizeof(double) * MIN(m,n));
  if (!S){

	LOG4CPLUS_ERROR(logger, "failed to allocate singular vectors.");
    return -2;
  }

  double * dummyU = NULL;

  // allocate array for VT
  int ldVT = MIN(m,n);
  double * VT = (double *) malloc (sizeof(double) * ldVT * n);
  if (!VT)
  {
	LOG4CPLUS_ERROR(logger, "failed to allocate VT.");
    return -2;
  }

  #ifdef __APPLE__
    #define DGESVD dgesvd_
    #define INTEGER long int
  #else
    #define DGESVD dgesvd_
    #define INTEGER long int
  #endif

  INTEGER M = m;
  INTEGER N = n;
  INTEGER LDA = ldA;
  INTEGER LDU = ldU;
  INTEGER LDVT = ldVT;
  INTEGER LWORK = lwork;
  INTEGER INFO = 0;

  LOG4CPLUS_TRACE(logger, "Calling LAPACK dgesvd routine...");
  //SUBROUTINE SGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT, WORK, LWORK, INFO )

  DGESVD (&jobu, &jobvt, &M, &N, A, &LDA,
	  S, dummyU, &LDU, VT,
	  &LDVT, work, &LWORK, &INFO);

  if (INFO != 0){

    int code = INFO;
	LOG4CPLUS_ERROR(logger, "Error: SVD solver returned non-zero exit code: " << code);
    free(VT);
    free(S);
    free(work);
    return code;
  }

  free(work);

  if (transpose)
  {
    InPlaceTransposeMatrix(m,n,VT);
    memcpy(A, VT, sizeof(double) * m * n);

    // swap m and n
    int bufferii = m;
    m = n;
    n = bufferii;
  }

  free(VT);

  //double totalEnergy = 0;
  std::ostringstream stringStream;
  for (int i=0; i< MIN(m,n); i++){
	stringStream << S[i] << " ";
  }
  LOG4CPLUS_DEBUG(logger, "Singular values: " << stringStream.str());

  // discard unneccesary modes
  LOG4CPLUS_TRACE(logger, "Discarding unnecessary components...");
  if (thresholdingSpecification->tresholdingType == ThresholdingSpecification::epsilonBased){
    DoTresholding_Epsilon(S,MIN(m,n),r,thresholdingSpecification->epsilon);
  }else{
    DoTresholding_NumberOfModes(S,MIN(m,n),r,thresholdingSpecification->rDesired);
  }
  // now, variable r has been set to the number of retained modes
  free(S);
  return 0;
}

