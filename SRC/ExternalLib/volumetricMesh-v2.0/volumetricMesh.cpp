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

#include <Log.h>
#include "volumetricMesh.h"
#include "abqparser.h"
#include "matrixMacros.h"
using namespace std;

// parses in the mesh, and returns the string corresponding to the element type
VolumetricMesh::VolumetricMesh(const char * filename, int numElementVertices, char ** elementType): numElementVertices_(numElementVertices)
{
  DEBUG_LOG("Parsing " << filename);

  // create buffer for element vertices
  int * v = (int*) malloc (sizeof(int) * numElementVertices_);

  // parse the abaqus file
  ABQParser abqParser;
  if (abqParser.open(filename) != 0){
	ERROR_LOG("could not open file: " << filename);
    free(v);
    throw 1;
  }

  int mode = 0;

  // pass 1: determine number of vertices and number of elements

  numElements_ = 0;
  numVertices_ = 0;

  numMaterials_ = 0;
  numElementSets_ = 0;
  numSolidSections_ = 0;
  *elementType = NULL;

  char s[4096];

  while (abqParser.getNextLine(s) != 0)
  {
    // now, the next input line is in s
    //printf("%s\n", s);

    // seek for *NODE
    if ((mode == 0) && (strncmp(s,"*NODE",5)==0))
    {
      mode = 1;
      continue;
    }

    // seek for *ELEMENT
    if ((mode == 1) && (strncmp(s,"*ELEMENT",8)==0))
    {
      mode = 2;
      // seek for "," and terminate elementType string there
      char * ch = &s[14];
      while ((*ch != ',') && (*ch != 0) && (*ch != ' '))
        ch++;
      *ch = 0;
      *elementType = (char*) malloc (sizeof(char) * (strlen(s) + 1 - 14));
      memcpy(*elementType, &s[14], sizeof(char) * (strlen(s) + 1 - 14));
      //*elementType = strdup(&s[14]);
      //printf("et: %p\n", *elementType);
      continue;
    }

    if ((mode == 2) && (s[0] == '*'))
    {
      mode = 3; // reached *ELSET
    }

    // in mode 1, increase numVertices counter
    if (mode == 1)
      numVertices_++;

    // in mode 2, read element info
    if (mode == 2)
      numElements_++;

    if (strncmp(s,"*MATERIAL,NAME=",15)==0)
    {
      numMaterials_++;
      continue;
    }

    if (strncmp(s,"*ELSET,ELSET=",13)==0)
    {
      numElementSets_++;
      continue;
    }

    if (strncmp(s,"*SOLIDSECTION,",14)==0)
    {
      numSolidSections_++;
      continue;
    }
  }

  abqParser.rewindToStart();

  if (numMaterials_ == 0)
  {
	WARN_LOG("no materials encountered in the input abaqus file.");
  }

  if (numSolidSections_ == 0)
  {
	WARN_LOG("no solid sections encountered in the input abaqus file.");
  }

  vertices_ = (Vec3d**) malloc (sizeof(Vec3d*) * numVertices_);
  elements_ = (MeshElement**) malloc (sizeof(MeshElement*) * numElements_);
  materials_ = (Material**) malloc (sizeof(Material*) * numMaterials_);
  elementSets_ = (ElementSet**) malloc (sizeof(ElementSet*) * numElementSets_);
  solidSections_ = (SolidSection**) malloc (sizeof(SolidSection*) * numSolidSections_);

  // pass 2: read in material data and element set data
  int numMaterials = -1;
  int numElementSets = -1;
  int numSolidSections = -1;
  mode = 0;
  char materialName[96];
  double E,nu,density;
  while (abqParser.getNextLine(s) != 0)
  {
    //printf("%s\n", s);
    int flag = 0;

    if ((mode == 11) && (s[0]=='*'))
      mode = 0;

    if ((mode == 0) && (strncmp(s,"*MATERIAL,NAME=",15)==0))
    {
      numMaterials++;
      strcpy(materialName,&s[15]);
      mode = 1;
      E = -1;
      nu = -1;
      density = -1;
      flag = 1;
    }

    if ((mode == 1) && (strncmp(s,"*ELASTIC",8)==0))
    {
      abqParser.getNextLine(s);
      sscanf(s,"%lf,%lf\n",&E,&nu);
      flag = 1;
    }

    if((mode == 1) && (strncmp(s,"*DENSITY",8)==0))
    {
      abqParser.getNextLine(s);
      sscanf(s,"%lf\n",&density);
      flag = 1;
    }

    if ((mode == 1) && (flag == 0))
    {
	  WARN_LOG("unknown material description encountered");
      mode = 0;
    }

    if ((mode == 1) && (E > 0) && (nu > 0) && (density > 0))
    {
      // create new material
      materials_[numMaterials] = new Material(E,nu,density,materialName);

      mode = 0;
      E = -1;
      nu = -1;
      density = -1;
    }

    if (mode == 12)
    {
      int start,end;
      sscanf(s,"%d,%d",&start,&end);
      for (int i=start; i<=end; i++)
        elementSets_[numElementSets]->addElement(i);
      mode=0;
    }

    if (mode == 11)
    {
      // we know that s[0] != '*' (i.e., not end of the list) as that case was already previously handled

      // parse the comma separated line
      char * pch;
      pch = strtok(s,",");
      while ((pch != NULL) && (isdigit(*pch)))
      {
        int newElement = atoi(pch);
        elementSets_[numElementSets]->addElement(newElement);
        pch = strtok (NULL, ",");
      }
    }

    if ((mode==0) && (strncmp(s,"*ELSET,ELSET=",13)==0))
    {
      numElementSets++;
      char elementSetName[96];
      strcpy(elementSetName,&s[13]);

      // check if elementSetName contains 'GENERATE'
      if (strstr(elementSetName,"GENERATE") == NULL)
      {
        // 'GENERATE' not found, read consecutive list to obtain member elements
        elementSets_[numElementSets] = new ElementSet(elementSetName);
        mode = 11;
      }
      else
      {
        // 'GENERATE' found, consecutive line will describe the set
        char * pos = strchr(elementSetName,','); // find comma before GENERATE
        *pos = '\0';
        elementSets_[numElementSets] = new ElementSet(elementSetName);
        mode = 12;
      }
    }

    if ((mode==0) && (strncmp(s,"*SOLIDSECTION,",14)==0))
    {
      numSolidSections++;
      char materialName[96];

      char * matPos = strstr(s,"MATERIAL=");
      strcpy(materialName,matPos+9);
      char * commaPos = strchr(materialName,',');
      *commaPos = '\0';

      char elementSetName[96];
      char * elementSetPos = strstr(s,"ELSET=");
      strcpy(elementSetName,elementSetPos+6);

      // seek for materialName
      int materialNum = -1;
      for(int material=0; material < numMaterials_; material++)
      {
        char * currentName = materials_[material]->name();
        if (strcmp(currentName,materialName) == 0)
        {
          materialNum = material;
          break;
        }
      }
      if (materialNum == -1)
      {
        printf("Error: material name %s not found among the materials.\n",materialName);
      }

      // seek for elementSetName
      int elementSetNum = -1;
      for(int elementSet=0; elementSet < numElementSets_; elementSet++)
      {
        char * currentElementSetName = elementSets_[elementSet]->name();
        if (strcmp(currentElementSetName,elementSetName) == 0)
        {
          elementSetNum = elementSet;
          break;
        }
      }
      if (elementSetNum == -1)
      {
        printf("Error: element set name %s not found among the element sets.\n",elementSetName);
      }

      solidSections_[numSolidSections] = new SolidSection(materialNum,elementSetNum);
    }
  }

  // pass 3: read vertices and elements

  abqParser.rewindToStart();

  numElements_ = 0;
  numVertices_ = 0;
  mode = 0;

  int numUnassignedElements = 0;
  while (abqParser.getNextLine(s) != 0)
  {
    //printf("%s\n", s);
    // now, next input line is in s

    // seek for *NODE
    if ((mode == 0) && (strncmp(s,"*NODE",5)==0))
    {
      mode = 1;
      continue;
    }

    // seek for *ELEMENT
    if ((mode == 1) && (strncmp(s,"*ELEMENT",8)==0))
    {
      mode = 2;
      continue;
    }

    if ((mode == 2) && (s[0] == '*'))
    {
      mode = 3; // reached *ELSET
    }

    // in mode 1, read the vertex info
    if (mode == 1)
    {
      int index;
      double x,y,z;
      sscanf(s,"%d,%lf,%lf,%lf",&index,&x,&y,&z);
      vertices_[numVertices_] = new Vec3d(x,y,z);

      numVertices_++;
    }

    // in mode 2, read element info
    if (mode == 2)
    {
      int index;
      sscanf(s,"%d",&index);
      char * ch = s;
      for(int i=0; i<numElementVertices_; i++)
      {
        // seek for next comma
        while((*ch != ',') && (*ch != 0))
          ch++;

        if (*ch == 0)
        {
		  ERROR_LOG("Error parsing line " << s);
          free(v);
          throw 12;
        }

        ch++;
        sscanf(ch,"%d",&v[i]);
      }
      
      for (int k=0; k<numElementVertices_; k++)
      {
        v[k]--; // compensate for fact that vertices are 1-numbered in abq files
      }

      // seek the element in all the solid sections
      int found = 0;
      int materialIndex = 0;
      for(int section=0; section < numSolidSections_; section++)
      {
        int elementSet = (solidSections_[section])->getElementSetIndex();
        int material = (solidSections_[section])->getMaterialIndex();

        // seek for element in elementSet
        if (elementSets_[elementSet]->isMember(index))
        {
          if (found == 0)
          {
            materialIndex = material;
            found = 1;
          }
          else
          {
			WARN_LOG("same element ("<<index<<") found in two solid sections.");
          }
        }
      }

      double E,nu,density;

      if (found == 0)
      {
        numUnassignedElements++;
        E = 1E4;
        nu = 0.45;
        density = 1000;
      }
      else
      {
        E = materials_[materialIndex]->E();
        nu = materials_[materialIndex]->nu();
        density = materials_[materialIndex]->density();
      }

      elements_[numElements_] = new MeshElement(E,nu,density,numElementVertices_,v);
      numElements_++;
    }
  }

  if (numUnassignedElements > 0){
	WARN_LOG(numUnassignedElements << " elements were not found in any of the solid sections");
	WARN_LOG("using default material parameters");
	WARN_LOG( "where E = "<<E <<" ,nu = "<<nu <<"and density = "<<density);
  }

  abqParser.close();

  free(v);
}

char * VolumetricMesh::getElementType(const char * filename, char ** elementType) 
{
  // parse the abaqus file
  ABQParser abqParser;
  *elementType = NULL;

  if (abqParser.open(filename) != 0){
    ERROR_LOG("could not open file: " << filename);
    return NULL;
  }

  char s[4096];
  while (abqParser.getNextLine(s) != 0)
  {
    // seek for *ELEMENT
    if (strncmp(s,"*ELEMENT",8)==0)
    {
      // seek for "," and terminate elementType string there
      char * ch = &s[14];
      while ((*ch != ',') && (*ch != 0) && (*ch != ' '))
        ch++;
      *ch = 0;
      *elementType = (char*) malloc (sizeof(char) * (strlen(s) + 1 - 14));
      memcpy(*elementType, &s[14], sizeof(char) * (strlen(s) + 1 - 14));
      break;
    }
  }

  abqParser.close();

  return *elementType;
}

VolumetricMesh::VolumetricMesh(int numVertices, double * vertices,
               int numElements, int numElementVertices, int * elements,
               double E, double nu, double density): numElementVertices_(numElementVertices)
{
  numElements_ = numElements;
  numVertices_ = numVertices;

  numMaterials_ = 1;
  numElementSets_ = 1;
  numSolidSections_ = 1;

  vertices_ = (Vec3d**) malloc (sizeof(Vec3d*) * numVertices_);
  elements_ = (MeshElement**) malloc (sizeof(MeshElement*) * numElements_);
  materials_ = (Material**) malloc (sizeof(Material*) * numMaterials_);
  elementSets_ = (ElementSet**) malloc (sizeof(ElementSet*) * numElementSets_);
  solidSections_ = (SolidSection**) malloc (sizeof(SolidSection*) * numSolidSections_);

  for(int i=0; i<numVertices; i++)
    vertices_[i] = new Vec3d(vertices[3*i+0], vertices[3*i+1], vertices[3*i+2]);

  Material * material = new Material(E, nu, density, "defaultMaterial");
  materials_[0] = material;

  ElementSet * elementSet = new ElementSet("allElements");

  int * v = (int*) malloc (sizeof(int) * numElementVertices_);
  for(int i=0; i<numElements; i++)
  {
    elementSet->addElement(i);
    for(int j=0; j<numElementVertices_; j++)
      v[j] = elements[numElementVertices_*i+j];
    elements_[i] = new MeshElement(E,nu,density,numElementVertices_,v);
  }
  free(v);

  elementSets_[0] = elementSet;

  SolidSection * solidSection = new SolidSection(0, 0);
  solidSections_[0] = solidSection;
}

VolumetricMesh::VolumetricMesh(VolumetricMesh & volumetricMesh, int numElements, int * elements, map<int,int> * vertexMap_)
{
  // determine vertices in the domain
  numElementVertices_ = volumetricMesh.numElementVertices();
  set<int> vertexSet;
  for(int i=0; i<numElements; i++)
    for(int j=0; j < numElementVertices_; j++)
      vertexSet.insert(volumetricMesh.vertexIndex(elements[i],j));

  // copy vertices into place and also into vertexMap
  numVertices_ = vertexSet.size();
  vertices_ = (Vec3d**) malloc (sizeof(Vec3d*) * numVertices_);
  set<int> :: iterator iter;
  int vertexNo = 0;
  map<int, int> vertexMap;
  for(iter = vertexSet.begin(); iter != vertexSet.end(); iter++)
  {
    vertices_[vertexNo] = new Vec3d(*(volumetricMesh.vertex(*iter)));
    vertexMap.insert(make_pair(*iter,vertexNo));
    vertexNo++;
  }

  if (vertexMap_ != NULL)
    *vertexMap_ = vertexMap;

  // copy elements
  numElements_ = numElements;
  elements_ = (MeshElement**) malloc (sizeof(MeshElement*) * numElements_);
  int * elementVertices = (int*) malloc (sizeof(int) * numElementVertices_);
  map<int,int> elementMap;
  for(int i=0; i < numElements; i++)
  {
    MeshElement * oldElement = volumetricMesh.element(elements[i]);
    for(int j=0; j< numElementVertices_; j++)
    {
      map<int,int> :: iterator iter2 = vertexMap.find(oldElement->vertex(j));
      if (iter2 == vertexMap.end())
      {
        printf("Internal error.\n");
        exit(1);
      }
      elementVertices[j] = iter2->second;
    }
    elements_[i] = new MeshElement(oldElement->E(), oldElement->nu(), oldElement->density(), numElementVertices_, elementVertices);
    elementMap.insert(make_pair(elements[i]+1, i+1)); // elements are 0-indexed
  }
  free(elementVertices);

  // copy materials
  numMaterials_ = volumetricMesh.numMaterials();
  numElementSets_ = volumetricMesh.numElementSets();
  numSolidSections_ = volumetricMesh.numSolidSections();

  materials_ = (Material**) malloc (sizeof(Material*) * numMaterials_);
  for(int i=0; i < numMaterials_; i++)
    materials_[i] = new Material(*(volumetricMesh.material(i)));

  // copy element sets; restrict element sets to the new mesh, also rename vertices to reflect new name
  vector<ElementSet*> elementSetsVector_;
  map<int,int> solidSectionsMap;
  for(int i=0; i < volumetricMesh.numElementSets(); i++)
  {
    ElementSet * oldElementSet = volumetricMesh.elementSet(i);
    set<int> oldElements;
    oldElementSet->getElements(oldElements);

    for(set<int> :: iterator iter = oldElements.begin(); iter != oldElements.end(); iter++)
    {
      if(*iter <= 0)
      {
        printf("Error 1.\n");
        exit(1);
      }
    }

    // construct the element list
    vector<int> newElements;
    for(set<int> :: iterator iter = oldElements.begin(); iter != oldElements.end(); iter++)
    {
      map<int,int> :: iterator iter2 = elementMap.find(*iter);
      if (iter2 != elementMap.end())
        newElements.push_back(iter2->second);
    }

    if (newElements.size() > 0)
    {
      ElementSet * newElementSet = new ElementSet(oldElementSet->name());
      for(unsigned int j=0; j<newElements.size(); j++)
      {
        if(newElements[j] <= 0)
        {
          printf("Error 2.\n");
          exit(1);
        }
        newElementSet->addElement(newElements[j]);
      }
      elementSetsVector_.push_back(newElementSet);
      solidSectionsMap.insert(make_pair(i,elementSetsVector_.size()-1));
    }
  }

  numElementSets_ = elementSetsVector_.size();
  elementSets_ = (ElementSet**) malloc (sizeof(ElementSet*) * numElementSets_);
  for(int j=0; j<numElementSets_; j++)
    elementSets_[j] = elementSets_[j];

  //printf("numElementSets_: %d\n",numElementSets_);

  // copy sections; remove empty ones
  vector<SolidSection*> sections_;
  for(int i=0; i < numSolidSections_; i++)
  {
    SolidSection * ssection = volumetricMesh.solidSection(i);
    map<int,int> :: iterator iter = solidSectionsMap.find(ssection->getElementSetIndex());
    if (iter != solidSectionsMap.end())
    {
      SolidSection * newSection = new SolidSection(ssection->getMaterialIndex(),iter->second);
      sections_.push_back(newSection);
    }
  }

  numSolidSections_ = sections_.size();
  solidSections_ = (SolidSection**) malloc (sizeof(SolidSection*) * numSolidSections_);
  for(int j=0; j<numSolidSections_; j++)
    solidSections_[j] = sections_[j];

  // sanity check
  // seek each element in all the solid sections
  for(int el=0; el<numElements_; el++)
  {
    int found = 0;
    for(int section=0; section < numSolidSections_; section++)
    {
      int elementSet = (solidSections_[section])->getElementSetIndex();

      // seek for element in elementSet
      if (elementSets_[elementSet]->isMember(el+1))
      {
        if (found != 0)
          printf("Warning: element %d is in more than one solid section.\n",el+1);
        else
          found = 1;
      }
    }
    if (found == 0)
      printf("Warning: element %d is not in any of the solid sections.\n",el+1);
  }

  // sanity check: make sure all elements are between bounds
  for(int i=0; i < numElementSets_; i++)
  {
    set<int> elts;
    elementSets_[i]->getElements(elts);
    for(set<int> :: iterator iter = elts.begin(); iter != elts.end(); iter++)
    {
      if (*iter <= 0)
        printf("Warning: encountered negative element index in element set %d.\n",i);
      if (*iter > numElements_)
        printf("Warning: encountered too large element index in element set %d.\n",i);
    }
  }
}

VolumetricMesh::~VolumetricMesh()
{
  for(int i=0; i< numVertices_; i++)
    delete(vertices_[i]);
  free(vertices_);

  for(int i=0; i< numElements_; i++)
    delete(elements_[i]);
  free(elements_);

  for(int i=0; i< numMaterials_; i++)
    delete(materials_[i]);
  free(materials_);
  
  for(int i=0; i< numElementSets_; i++)
    delete(elementSets_[i]);
  free(elementSets_);

  for(int i=0; i< numSolidSections_; i++)
    delete(solidSections_[i]);
  free(solidSections_);
}

int VolumetricMesh::save(char * filename, char * elementName) // saves the mesh to an abq file
{       
  FILE * fout = fopen(filename,"wa");
  if (!fout)
  {       
    printf("Error: could not write to %s.\n",filename);
    return 1;
  }         
          
  fprintf(fout,"*NODE\n");
          
  for(int i=0; i < numVertices_; i++)  {
    Vec3d v = *vertex(i);
    fprintf(fout,"%d, %.15f,  %.15f,  %.15f\n",i+1,v[0],v[1],v[2]);
  }   

  fprintf(fout,"*ELEMENT, TYPE=%s\n", elementName);
  for(int el=0; el < numElements_; el++)
  {   
    fprintf(fout,"%d, ",el+1);
    for(int j=0; j < numElementVertices_; j++)
    {   
      fprintf(fout,"%d",vertexIndex(el,j)+1);
      if (j != numElementVertices_ - 1)
        fprintf(fout,", ");
    } 
    fprintf(fout,"\n");
  }     
        
  fprintf(fout,"*ELSET,ELSET=EALL,GENERATE\n");
  fprintf(fout,"  1,%d\n",numElements_);

  fclose(fout);
  return 0;
}   

Vec3d VolumetricMesh::elementCenter(int el)
const{
  Vec3d pos(0,0,0);
  for(int i=0; i<numElementVertices_;i++)
    pos += *vertex(el,i);

  pos *= 1.0 / numElementVertices_;

  return pos;
}

void VolumetricMesh::getVerticesInElements(vector<int> & elements, vector<int> & vertices)
{
  set<int> ver;
  for(unsigned int i=0; i< elements.size(); i++)
    for(int j=0; j< numElementVertices_; j++)
      ver.insert(vertexIndex(elements[i],j));

  vertices.clear();
  set<int>::iterator iter;
  for(iter = ver.begin(); iter != ver.end(); iter++)
    vertices.push_back(*iter);
}

void VolumetricMesh::getElementsTouchingVertices(vector<int> & vertices, vector<int> & elements)
{
  set<int> ver;
  for(unsigned int i=0; i<vertices.size(); i++)
    ver.insert(vertices[i]);

  elements.clear();
  for(int i=0; i< numElements_; i++)
  {
    set<int> :: iterator iter;
    for(int j=0; j<numElementVertices_; j++)
    {
      iter = ver.find(vertexIndex(i,j));
      if (iter != ver.end())
      {
        elements.push_back(i);
        break;
      }
    }
  }
}

void VolumetricMesh::getVertexNeighborhood(vector<int> & vertices, vector<int> & neighborhood)
{
  vector<int> elements;
  getElementsTouchingVertices(vertices,elements);
  getVerticesInElements(elements,neighborhood);
}

void VolumetricMesh::getInertiaParameters(double & mass, Vec3d & centerOfMass, Mat3d & inertiaTensor)
{
  mass = 0.0;
  centerOfMass[0] = centerOfMass[1] = centerOfMass[2] = 0;
  inertiaTensor[0][0] = inertiaTensor[0][1] = inertiaTensor[0][2] = 0;
  inertiaTensor[1][0] = inertiaTensor[1][1] = inertiaTensor[1][2] = 0;
  inertiaTensor[2][0] = inertiaTensor[2][1] = inertiaTensor[2][2] = 0;

  // compute mass, center of mass, inertia tensor
  for(int i=0; i< numSolidSections(); i++)
  {
    SolidSection * solidSection_ = solidSection(i);
    ElementSet * elementSet_ = elementSet(solidSection_->getElementSetIndex());
    double density_ = material(solidSection_->getMaterialIndex())->density();

    set<int> elements;
    elementSet_->getElements(elements);

    // over all elements in the section
    for(set<int> :: iterator iter = elements.begin(); iter != elements.end(); iter++)
    {
      int element = *iter - 1;
      double elementVolume = getElementVolume(element);
      double elementMass = elementVolume * density_;
      //printf("elementMass: %G section size: %d\n",elementMass,elements.size());

      mass += elementMass;
      Vec3d elementCenter_ = elementCenter(element);
      centerOfMass += elementMass * elementCenter_;

      Mat3d elementITUnitDensity;
      getElementInertiaTensor(element, elementITUnitDensity);

      double a = elementCenter_[0];
      double b = elementCenter_[1];
      double c = elementCenter_[2];

      Mat3d elementITCorrection
       (  b*b + c*c, -a*b, -a*c,
         -a*b, a*a + c*c, -b*c,
         -a*c, -b*c, a*a + b*b );

      Mat3d elementIT = density_ * elementITUnitDensity + elementMass * elementITCorrection;

      inertiaTensor += elementIT;
    }
  }

  //printf("final mass: %G\n",mass);
  centerOfMass /= mass;

  // correct inertia tensor so it's around the center of mass
  double a = centerOfMass[0];
  double b = centerOfMass[1];
  double c = centerOfMass[2];

  Mat3d correction
       ( b*b + c*c, -a*b, -a*c,
         -a*b, a*a + c*c, -b*c,
         -a*c, -b*c, a*a + b*b );

  inertiaTensor -= mass * correction;
}

void VolumetricMesh::getMeshGeometricParameters(double * centerX, double * centerY, double * centerZ, double * radius)
{
  *centerX = 0.0;
  *centerY = 0.0;
  *centerZ = 0.0;

  for(int i=0; i < numVertices_ ; i++)
  {
    Vec3d * vertex_ = vertex(i);
    *centerX += (*vertex_)[0];
    *centerY += (*vertex_)[1];
    *centerZ += (*vertex_)[2];
  }

  *centerX /= numVertices_;
  *centerY /= numVertices_;
  *centerZ /= numVertices_;

  *radius = 0;
  Vec3d center(*centerX,*centerY,*centerZ);

  for(int i=0; i < numVertices_ ; i++)
  {
    Vec3d * vertex_ = vertex(i);

    double dist = len(*vertex_-center);

    if (dist > *radius)
      *radius = dist;
  }
}

int VolumetricMesh::getClosestVertex(Vec3d pos)
{
  // linear scan
  double closestDist = DBL_MAX;
  int closestVertex = -1;

  for(int i=0; i<numVertices_; i++)
  {
    Vec3d * vertexPosition = vertices_[i];
    double dist = len(pos - *vertexPosition);
    if (dist < closestDist)
    {
      closestDist = dist;
      closestVertex = i;
    }
  }

  return closestVertex;
}

int VolumetricMesh::getClosestElement(Vec3d pos)
{
  // linear scan
  double closestDist = DBL_MAX;
  int closestElement = 0;
  for(int element=0; element < numElements_; element++)
  {
    Vec3d center = elementCenter(element);
    double dist = len(pos - center);
    if (dist < closestDist)
    {
      closestDist = dist;
      closestElement = element;
    }
  }

  return closestElement;
}

int VolumetricMesh::getContainingElement(Vec3d pos)
{
  // linear scan
  for(int element=0; element < numElements_; element++)
  {
    if (containsVertex(element, pos))
      return element;
  }

  return -1;
}

void VolumetricMesh::setSingleMaterial(double E, double nu, double density)
{
  // erase previous materials
  for(int i=0; i<numMaterials_; i++)
    free(materials_[i]);
  free(materials_);

  for(int i=0; i<numElementSets_; i++)
    free(elementSets_[i]);
  free(elementSets_);

  for(int i=0; i<numSolidSections_; i++)
    free(solidSections_[i]);
  free(solidSections_);

  // add a single material
  numMaterials_ = 1;
  numElementSets_ = 1;
  numSolidSections_ = 1;

  materials_ = (Material**) malloc (sizeof(Material*) * numMaterials_);
  elementSets_ = (ElementSet**) malloc (sizeof(ElementSet*) * numElementSets_);
  solidSections_ = (SolidSection**) malloc (sizeof(SolidSection*) * numSolidSections_);

  Material * material = new Material(E, nu, density, "defaultMaterial");
  materials_[0] = material;

  ElementSet * elementSet = new ElementSet("allElements");

  int * v = (int*) malloc (sizeof(int) * numElementVertices_);
  for(int i=0; i<numElements_; i++)
  {
    for(int j=0; j<numElementVertices_; j++)
      v[j] = vertexIndex(i,j);
    free(elements_[i]);
    elements_[i] = new MeshElement(E,nu,density,numElementVertices_,v);
    elementSet->addElement(i+1);
  }
  free(v);

  elementSets_[0] = elementSet;

  SolidSection * solidSection = new SolidSection(0, 0);
  solidSections_[0] = solidSection;
}

void VolumetricMesh::PropagateMaterialsToElements()
{
  for(int sectionIndex=0; sectionIndex < numSolidSections_; sectionIndex++)
  {
    SolidSection * section = solidSections_[sectionIndex];

    int materialIndex = section->getMaterialIndex();
    double E = materials_[materialIndex]->E();
    double nu = materials_[materialIndex]->nu();
    double density = materials_[materialIndex]->density();

    int elementSetIndex = section->getElementSetIndex();
    set<int> setElements;
    elementSets_[elementSetIndex]->getElements(setElements);

    int eltCounter = 0;
    for(set<int> :: iterator iter = setElements.begin(); iter != setElements.end(); iter++)
    {
      //printf("elt=%d %p\n", *iter, elements_[*iter]); fflush(NULL);
      int elt = *iter - 1;
      elements_[elt]->setE(E);
      elements_[elt]->setNu(nu);
      elements_[elt]->setDensity(density);
      eltCounter++;
    }
  }
}

int VolumetricMesh::buildInterpolationWeights(int numTargetLocations, 
	const float * targetLocations, int ** vertices, double ** weights, 
	double zeroThreshold, vector<int> * elementList, int verbose){

  // modified by lsw
  TRACE_FUN();

  // allocate interpolation arrays  
  *vertices = new int[numElementVertices_ * numTargetLocations];
  *weights = new double[numElementVertices_ * numTargetLocations];
  double * baryWeights = (double*) malloc (sizeof(double) * numElementVertices_);
  int numExternalVertices = 0;

  INFO_LOG("begin to build the interpolate weights, which needs a long time please wait ...");
  INFO_LOG("it will stop when counter reach " << numTargetLocations);
  for (int i=0; i < numTargetLocations; i++){// over all interpolation locations

    if ((verbose) && (i % 100 == 0)){
      printf("%d ", i); fflush(NULL);
    }
    const Vec3d pos = Vec3d(targetLocations[3*i+0],targetLocations[3*i+1],targetLocations[3*i+2]);

    // find element containing pos
    int element = getContainingElement(pos);

    if (element < 0){

      element = getClosestElement(pos);
      numExternalVertices++;
    }

    if (elementList != NULL){

      elementList->push_back(element);
	}

    computeBarycentricWeights(element, pos, baryWeights);

    if (zeroThreshold > 0){
      // check whether vertex is close enough to the sqmesh
      double minDistance = DBL_MAX;
      int numElementVertices_ = numElementVertices();
      int assignedZero = 0;
      for(int ii=0; ii< numElementVertices_; ii++){
        Vec3d * vpos = vertex(element, ii);
        if (len(*vpos-pos) < minDistance){
          minDistance = len(*vpos-pos);
        }
      }

      if (minDistance > zeroThreshold){
        // assign zero weights
        for(int ii=0; ii < numElementVertices_; ii++){
          baryWeights[ii] = 0.0;
		}
        assignedZero++;
        continue;
      }
    }

    for(int ii=0; ii<numElementVertices_; ii++){
      (*vertices)[numElementVertices_ * i + ii] = vertexIndex(element, ii);
      (*weights)[numElementVertices_ * i + ii] = baryWeights[ii];
    }
  }
  free(baryWeights);
  return numExternalVertices;
}

int VolumetricMesh::getNumInterpolationElementVertices(char * filename){

  FILE * fin = fopen(filename, "ra");
  if (!fin)
  {
    printf("Error: unable to open file %s.\n", filename);
    return -1;
  }

  char s[4096];
  if (fgets(s, 4096, fin) == NULL)
  {
    printf("Error: incorrect first line of file %s.\n", filename);
    return -2;
  }
  fclose(fin);

  ABQParser::beautifyLine(s, 1);

  int slen = strlen(s);
  int count = 0;
  for(int i=0; i<slen; i++)
    if (s[i] == ' ')
      count++;

  if (count % 2 == 1)
  {
    printf("Error: odd number of whitespaces in the first line of file %s.\n", filename);
    return -3;
  }

  return count / 2;
}

int VolumetricMesh::loadInterpolationWeights(const char * filename, int numTargetLocations, 
											 int numElementVertices_, int ** vertices, double ** weights){
  FILE * fin = fopen(filename, "ra");
  if (!fin){
    printf("Error: unable to open file %s.\n", filename);
    return 2;
  }

  // allocate interpolation arrays
  // modified by lsw
  *vertices = new int[numElementVertices_ * numTargetLocations];
  *weights = new double[numElementVertices_ * numTargetLocations];

  int numVertices = -1;
  int currentVertex;

  // read the elements one by one and accumulate entries
  while (numVertices < numTargetLocations-1)
  {
    numVertices++;

    if (feof(fin))
    {
      printf("Error: interpolation file is too short. Num vertices in interp file: %d . Should be: %d .\n", numVertices, numTargetLocations);
      return 1;
    }

    fscanf(fin, "%d", &currentVertex);
    if (currentVertex != numVertices)
    {
      printf("Error: consecutive vertex index at position %d mismatch.\n", currentVertex);
      return 1;
    }

    for(int j=0; j<numElementVertices_; j++)
      fscanf(fin,"%d %lf", &((*vertices)[currentVertex * numElementVertices_ + j]), &((*weights)[currentVertex * numElementVertices_ + j]) );

    fscanf(fin,"\n");
  }

  fclose(fin);

  return 0;
}

int VolumetricMesh::saveInterpolationWeights(const char * filename, int numTargetLocations, 
											 int numElementVertices_, const int * vertices, const double * weights){

  FILE * fin = fopen(filename, "wa");
  if (!fin){
	
    printf("Error: unable to open file %s.\n", filename);
    return 1;
  }

  // read the elements one by one and accumulate entries
  for(int currentVertex=0; currentVertex < numTargetLocations; currentVertex++){

    fprintf(fin, "%d", currentVertex);
    for(int j=0; j<numElementVertices_; j++)
      fprintf(fin," %d %lf", 
        vertices[currentVertex * numElementVertices_ + j], 
        weights[currentVertex * numElementVertices_ + j]);

    fprintf(fin,"\n");
  }

  fclose(fin);

  return 0;
}

void VolumetricMesh::interpolate(const double * u, float * uTarget, 
								 int numTargetLocations, int numElementVertices_, 
								 const int * vertices, const double * weights){

  for(int i=0; i< numTargetLocations; i++){

    Vec3d defo(0,0,0);
    for(int j=0; j<numElementVertices_; j++){
      const int k = vertices[numElementVertices_ * i + j];// volumetricMeshVertexIndex.
      const Vec3d volumetricMeshVertexDefo = Vec3d(u[3*k+0],u[3*k+1],u[3*k+2]);
      defo += weights[numElementVertices_ * i + j] * volumetricMeshVertexDefo;
    }
    uTarget[3*i+0] = defo[0];
    uTarget[3*i+1] = defo[1];
    uTarget[3*i+2] = defo[2];
  }
}

void VolumetricMesh::dumpMeshData(int * numVertices, double ** vertices, int * numElements, 
								  int * numElementVertices, int ** elements){
  *numVertices = numVertices_;
  *numElements = numElements_;
  *numElementVertices = numElementVertices_;

  *vertices = (double*) malloc (sizeof(double) * 3 * numVertices_);
  *elements = (int*) malloc (sizeof(int) * numElementVertices_ * numElements_);

  for(int i=0; i<numVertices_; i++)
  {
    Vec3d v = *vertex(i);
    (*vertices)[3*i+0] = v[0];
    (*vertices)[3*i+1] = v[1];
    (*vertices)[3*i+2] = v[2];
  }

  for(int i=0; i<numElements_; i++)
  {
    MeshElement * elt = element(i);
    for(int j=0; j<numElementVertices_; j++)
      (*elements)[numElementVertices_ * i + j] = elt->vertex(j);
  }
}

void VolumetricMesh::computeGravity(double * gravityForce, double g, bool addForce){
  if (!addForce)
    memset(gravityForce, 0, sizeof(double) * 3 * numVertices_);

  double invNumElementVertices = 1.0 / numElementVertices();

  for(int el=0; el < numElements_; el++)
  {
    double volume = getElementVolume(el);
    double density = element(el)->density();
    double mass = density * volume;
    for(int j=0; j<numElementVertices(); j++)
      gravityForce[3 * (vertexIndex(el,j)) + 1] -= invNumElementVertices * mass * g;
  }  
}

void VolumetricMesh::applyDeformation(double * u){
  for(int i=0; i<numVertices_; i++)
  {
    Vec3d * v = vertex(i);
    (*v)[0] += u[3*i+0];
    (*v)[1] += u[3*i+1];
    (*v)[2] += u[3*i+2];
  }
}

void VolumetricMesh::subsetMesh(std::set<int> & elements, int removeIsolatedVertices){
  int numRemovedElements = 0;
  for(int el=0; el<numElements_; el++)
  {
    if (elements.find(el) == elements.end())
    {
      delete(elements_[el]);
      elements_[el] = NULL;
      numRemovedElements++;
    }
  }

  int head = 0;
  int tail = 0;

  while (tail < numElements_)
  {
    if (elements_[tail] != NULL)
    {
      elements_[head] = elements_[tail];
      head++;
    }
    tail++;
  }
   
  numElements_ -= numRemovedElements;
  elements_ = (MeshElement**) realloc (elements_, sizeof(MeshElement*) * numElements_ );

  if (removeIsolatedVertices)
  {
    set<int> retainedVertices;
    for(int el=0; el<numElements_; el++)
      for(int j=0; j < numElementVertices_; j++)
        retainedVertices.insert(vertexIndex(el,j));

    int numRemovedVertices = 0;
    for(int v=0; v<numVertices_; v++)
    {
      if (retainedVertices.find(v) == retainedVertices.end())
      {
        delete(vertices_[v]);
        vertices_[v] = NULL;
        numRemovedVertices++;
      }
    }
  
    int head = 0;
    int tail = 0;
  
    int * renamingFunction = (int*) malloc (sizeof(int) * numVertices_);
    while (tail < numVertices_)
    {
      if (vertices_[tail] != NULL)
      {
        renamingFunction[tail] = head;
        vertices_[head] = vertices_[tail];
        head++;
      }
      tail++;
    }

    // rename vertices inside the elements
    for(int el=0; el<numElements_; el++)
      for(int j=0; j < numElementVertices_; j++)
        elements_[el]->setVertex(j, renamingFunction[vertexIndex(el,j)]);
   
    free(renamingFunction);
    numVertices_ -= numRemovedVertices;
    vertices_ = (Vec3d**) realloc (vertices_, sizeof(Vec3d*) * numVertices_ );
  }
}

int VolumetricMesh::exportToEle(char * baseFilename, int includeSolidSections){
  char s[4096];
  sprintf(s, "%s.ele", baseFilename);

  FILE * fout = fopen(s,"wa");
  if (!fout)
  {       
    printf("Error: could not write to %s.\n",s);
    return 1;
  }         

  int * elementSection = NULL;
  if (includeSolidSections)
  {
    elementSection = (int*) malloc (sizeof(int) * numElements());
    for(int el=0; el<numElements(); el++)
    {
      int found = 0;
      for(int section=0; section < numSolidSections_; section++)
      {
        int elementSet = (solidSections_[section])->getElementSetIndex();

        // seek for element in elementSet
        if (elementSets_[elementSet]->isMember(el+1))
        {
          if (found != 0)
            printf("Warning: element %d is in more than one solid section.\n",el+1);
          else
            found = section+1;
        }
      }
      if (found == 0)
        printf("Warning: element %d is not in any of the solid sections.\n",el+1);
      elementSection[el] = found;
    }
  }

  if (includeSolidSections)
    fprintf(fout,"%d %d %d\n", numElements_, numElementVertices_, 1);
  else
    fprintf(fout,"%d %d %d\n", numElements_, numElementVertices_, 0);

  for(int el=0; el < numElements_; el++)
  {   
    fprintf(fout,"%d ",el+1);
    for(int j=0; j < numElementVertices_; j++)
    {   
      fprintf(fout,"%d",vertexIndex(el,j)+1);
      if (j != numElementVertices_ - 1)
        fprintf(fout," ");
    } 
    if (includeSolidSections)
    {
      fprintf(fout," %d", elementSection[el]);
    }
    fprintf(fout,"\n");
  }     
        
  fprintf(fout,"# generated by the volumetricMesh class\n");

  fclose(fout);

  if (includeSolidSections)
    free(elementSection);

  sprintf(s, "%s.node", baseFilename);

  fout = fopen(s,"wa");
  if (!fout)
  {       
    printf("Error: could not write to %s.\n",s);
    return 1;
  }         
          
  fprintf(fout,"%d %d %d %d\n", numVertices_, 3, 0, 0);
  for(int v=0; v < numVertices_; v++)
  {   
    fprintf(fout,"%d ",v+1);
    Vec3d vtx = *vertex(v);
    fprintf(fout,"%.15f %.15f %.15f\n", vtx[0], vtx[1], vtx[2]);
  }     
        
  fprintf(fout,"# generated by the volumetricMesh class\n");

  fclose(fout);

  return 0;
}

// added by lsw
void VolumetricMesh::vertices(std::vector<double> &v)const{
  
  assert (this->numVertices() >= 0);

  v.resize(this->numVertices()*3);
  for (int i = 0; i < this->numVertices(); ++i){

	assert (vertices_ != NULL);
    v[i*3+0] = (*vertices_[i])[0];
    v[i*3+1] = (*vertices_[i])[1];
    v[i*3+2] = (*vertices_[i])[2];
  }
}

void VolumetricMesh::setMaterial(int i, Material & material) { 

  *(materials_[i]) = material;
  printf("pp: %p\n", elements_[50]); 
  PropagateMaterialsToElements(); 
}

void VolumetricMesh::barycentersOfGroup(const std::vector<std::set<int> >&g, double *bary_p)const{

  assert (bary_p != NULL);
  for (int i = 0; i < (int)g.size(); ++i){

	bary_p[i*3+0] = 0;
	bary_p[i*3+1] = 0;
	bary_p[i*3+2] = 0;
	
	const set<int> &one_group = g[i];
	set<int>::const_iterator it = one_group.begin();
	for (; it != one_group.end(); ++it){

	  const Vec3d &v = *(this->vertex(*it));
	  bary_p[i*3+0] += v[0];
	  bary_p[i*3+1] += v[1];
	  bary_p[i*3+2] += v[2];
	}
	if (one_group.size() > 0){

	  bary_p[i*3+0] /= one_group.size();
	  bary_p[i*3+1] /= one_group.size();
	  bary_p[i*3+2] /= one_group.size();
	}
  }
}

void VolumetricMesh::barycentersOfGroup(const std::vector<int>&one_group, double bary_p[3])const{

  bary_p[0] = 0;
  bary_p[1] = 0;
  bary_p[2] = 0;
  
  for (int i = 0; i < (int)one_group.size(); ++i){

	const Vec3d &v = *(this->vertex(i));
	bary_p[0] += v[0];
	bary_p[1] += v[1];
	bary_p[2] += v[2];
  }
  if (one_group.size() > 0){

	bary_p[0] /= one_group.size();
	bary_p[1] /= one_group.size();
	bary_p[2] /= one_group.size();
  }
}

void VolumetricMesh::barycentersOfGroup(const std::set<int>&one_group, double bary_p[3])const{
  
  bary_p[0] = 0;
  bary_p[1] = 0;
  bary_p[2] = 0;
  
  std::set<int>::const_iterator it = one_group.begin();
  for (; it != one_group.end(); it++ ){

	const Vec3d &v = *(this->vertex(*it));
	bary_p[0] += v[0];
	bary_p[1] += v[1];
	bary_p[2] += v[2];
  }

  if (one_group.size() > 0){

	bary_p[0] /= one_group.size();
	bary_p[1] /= one_group.size();
	bary_p[2] /= one_group.size();
  }
}
