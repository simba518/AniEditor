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

#ifndef _ELEMENTSET_H_
#define _ELEMENTSET_H_

#include "string.h"
#include <set>

// a set of elements (but can store any set; e.g., vertices)

class ElementSet
{
public:

  ElementSet(char * name);
  ElementSet(ElementSet & elementSet);

  inline void addElement(int element);
  inline void clear();
  inline void getElements(std::set<int> & elements);
  inline bool isMember(int element);

  inline char * name();

protected:
  char name_[24];
  std::set<int> elements;
};

inline ElementSet::ElementSet(char * name) { strcpy(name_,name); }
inline ElementSet::ElementSet(ElementSet & elementSet) { elements = elementSet.elements; strcpy(name_, elementSet.name()); }
inline void ElementSet::addElement(int element) { elements.insert(element); }
inline void ElementSet::clear() { elements.clear(); }
inline void ElementSet::getElements(std::set<int> & elements) { elements = this->elements; }
inline bool ElementSet::isMember(int element) {return (elements.find(element) != elements.end());}
inline char * ElementSet::name() { return name_; }

#endif

