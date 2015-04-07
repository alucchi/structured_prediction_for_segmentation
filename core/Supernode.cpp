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

// standard libraries
#include <vector>

// SliceMe
#include "Supernode.h"

using namespace std;

//------------------------------------------------------------------------------

double node_square_distance(const node& a, const node& b) {
  double dist = 0;
  double t = (a.x-b.x);
  dist += t*t;
  t = (a.y-b.y);
  dist += t*t;
  t = (a.z-b.z);
  dist += t*t;
  return dist;
}

//------------------------------------------------------------------------------

uchar supernode::getLabel() {
  if(data) {
    return data->label;
  } else {
    return INVALID_LABEL;
  }
}

void supernode::setLabel(uchar _label)
{
  if(!data) {
    data = new supernode_data;
  }

  data->label = _label;
}

void supernode::setData(uchar _label, int _nr_class, probType* _prob_estimates)
{
  if(!data) {
    data = new supernode_data(_nr_class, _prob_estimates);
  }
  /*
    else {
    printf("[Supernode] Warning : data has already been set for supernode %d\n", id);
  }
  */

  data->label = _label;
}

/**
 * Get center of the supernode with corresponding given id
 * @param center will be initialized by this function
 */
void supernode::getCenter(node& center)
{
  //TODO : might want to check that center is not outside supernode
  ulong cx = 0;
  ulong cy = 0;
  ulong cz = 0;
  node n;
  nodeIterator ni = getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    cx += n.x;
    cy += n.y;
    cz += n.z;
    ni.next();
  }
  int _size = size();
  if(_size != 0) {
    center.x = cx/_size;
    center.y = cy/_size;
    center.z = cz/_size;
  }
}

uint supernode::size()
{
  uint nodeSize = 0;
  // count number of pixels in line containers
  for(vector<lineContainer*>::iterator it = lines.begin();
      it != lines.end(); it++)
    nodeSize += (*it)->length;

  return nodeSize + nodes.size();
}

