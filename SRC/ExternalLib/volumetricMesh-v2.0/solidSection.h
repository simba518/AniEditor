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

#ifndef _SOLIDSECTION_H_
#define _SOLIDSECTION_H_

// a volumetric mesh solid section
// i.e., a set of elements sharing the same material

class SolidSection
{
public:
  SolidSection(int materialIndex, int elementSetIndex);
  inline int getMaterialIndex();
  inline int getElementSetIndex();

protected:
  int materialIndex_, elementSetIndex_;
};

inline SolidSection::SolidSection(int materialIndex, int elementSetIndex):
    materialIndex_(materialIndex),elementSetIndex_(elementSetIndex) {}
inline int SolidSection::getMaterialIndex() { return materialIndex_; }
inline int SolidSection::getElementSetIndex() { return elementSetIndex_; }

#endif

