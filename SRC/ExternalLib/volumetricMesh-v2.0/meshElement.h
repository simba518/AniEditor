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

#ifndef _MESHELEMENT_H_
#define _MESHELEMENT_H_

// a single mesh element: vertices, and its material properties

class MeshElement
{
public:
  MeshElement(double E, double nu, double density, int numVertices, int * vertices);
  MeshElement(MeshElement & meshElement);
  ~MeshElement();

  // Young's modulus, Poisson ration and mass density, 
  // and Lame's constants lambda and mu
  inline double E() { return E_; }
  inline double nu() { return nu_; }
  inline double density() { return density_; }
  inline double lambda() { return (nu_ * E_) / (1+nu_) / (1-2*nu_); }
  inline double mu() { return E_ / (2*(1+nu_)); }

  inline int vertex(int i) { return vertices_[i]; }
  inline int vertex(int i)const { return vertices_[i]; } // add by lsw
  inline int * vertices() { return vertices_; }
  inline const int * vertices() const { return vertices_; }

  inline void setVertex(int i, int vertex) { vertices_[i] = vertex; }

  inline void setE(double E) { E_ = E; }
  inline void setNu(double nu) { nu_ = nu; }
  inline void setDensity(double density) { density_ = density; }

protected:
  double E_, nu_;
  double density_;
  int numVertices_;
  int * vertices_;
};

#endif

