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

#ifndef _MATERIAL_H_
#define _MATERIAL_H_

#include "string.h"

// stores a StVK material (E, nu, density)

class Material
{
public:
  Material(double E, double nu, double density, char * name): E_(E), nu_(nu), density_(density) { strcpy(name_,name); }
  Material(Material & material) { E_ = material.E(); nu_ = material.nu(); density_ = material.density(); strcpy(name_,material.name()); }

  inline double E() { return E_; } // Young's modulus
  inline double nu() { return nu_; }  // Poisson ratio
  inline double density() { return density_; } // density
  inline char * name() { return name_; }  // material name

  // warning: altering these alone will not propagate the information to 
  // elements inside VolumetricMesh 
  // (must call "VolumetricMesh::PropagateMaterialsToElements")
  inline void setE(double E) { E_ = E; }
  inline void setNu(double nu) { nu_ = nu; }
  inline void setDensity(double density) { density_ = density; }
  inline void setName(char * name) { strcpy(name_,name); }

protected:
  double E_,nu_,density_;
  char name_[24];
};

#endif

