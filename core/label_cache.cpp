/////////////////////////////////////////////////////////////////////////
// This program is free software; you can redistribute it and/or       //
// modify it under the terms of the GNU General Public License         //
// version 2 as published by the Free Software Foundation.             //
//                                                                     //
// This program is distributed in the hope that it will be useful, but //
// WITHOUT ANY WARRANTY; without even the implied warranty of          //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   //
// General Public License for more details.                            //
//                                                                     //
// Written and (C) by Aurelien Lucchi                                  //
// Contact <aurelien.lucchi@gmail.com> for comments & bug reports      //
/////////////////////////////////////////////////////////////////////////

#include "label_cache.h"

LabelCache* LabelCache::pInstance = 0; // initialize pointer

LabelCache::~LabelCache()
{
  clear();
}

void LabelCache::clear()
{
  for(map<int, LABEL>::iterator it = labels.begin();
      it != labels.end(); ++it) {
    delete[] it->second.nodeLabels;
  }
  labels.clear();
}

bool LabelCache::exists(int id)
{
  bool labelFound = false;
  map<int, LABEL>::iterator lookup = labels.find(id);
  if(lookup != labels.end()) {
    labelFound = true;
  }
  return labelFound;
}

bool LabelCache::getLabel(int id, LABEL& l)
{
  bool labelFound = false;
  map<int, LABEL>::iterator lookup = labels.find(id);
  if(lookup != labels.end()) {
    l = labels[id];
    labelFound = true;
  }
  return labelFound;
}

void LabelCache::setLabel(int id, LABEL& l)
{
  labels[id] = l;
}
