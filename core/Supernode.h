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

#ifndef SUPERNODE_H
#define SUPERNODE_H

// standard libraries
#include <stdio.h>
#include <string.h>

// SliceMe
#include "globalsE.h"

//------------------------------------------------------------------------------

//typedef unsigned short sidType;
typedef int sidType; //should be signed (some functions like getSupernodeLabel return -1)
typedef uchar labelType;
typedef float probType;

//------------------------------------------------------------------------------

/**
 * Node structure using bit-fields to reduce amount of memory
 * A node is stored on 64 bits (could be reduced to 32 bits for 32bits machines)
 */
struct node
{
  int x:22; // maximum width  = 4194304
  int y:22; // maximum height = 4194304
  int z:20; // maximum depth  = 1048576

  node(){ x = 0; y = 0; z = 0; }
};

struct lineContainer
{
  node coord;
  uint length;

  lineContainer() {
    length = 0;
  }
};

//------------------------------------------------------------------------------

class supernode_data
{
 public:
  uchar label; // background/foreground...
  probType* prob_estimates;

  /**
   * @param _prob_size is the size of the probability vector (=number of classes in general)
   */
  supernode_data() {
    prob_estimates = 0;
  }

  /**
   * @param _prob_size is the size of the probability vector (=number of classes in general)
   */
  supernode_data(int _prob_size, probType* _prob_estimates) {
    prob_estimates = new probType[_prob_size];
    if(_prob_estimates)
      memcpy(prob_estimates,_prob_estimates,_prob_size*sizeof(probType));
    else
      memset(prob_estimates,0,_prob_size*sizeof(probType));
  }

  ~supernode_data() {
    delete[] prob_estimates;
  }
};


/**
 * Class to iterate over nodes
 * 2 types of storage are available : nodes and lines (i.e run length encoding)
 */
class nodeIterator
{
 public:
  nodeIterator(vector<lineContainer*>* _lines,
               vector<node*>* _nodes)
    {
      lines = _lines;
      nodes = _nodes;
      goToBegin();
    }

  ~nodeIterator()
    {
    }

  void goToBegin()
  {
    itLines = lines->begin();
    nodeIdx = 0;
    lineIdx = 1;
    itNodes = nodes->begin();
  }

  bool isAtEnd()
  {
    if(itLines == lines->end() && itNodes == nodes->end())
      return true;
    else
      return false;
  }

  inline void next()
  {
    if(itLines != lines->end())
      {
        //printf("next (*itLines)->length %d %d\n", nodeIdx, (*itLines)->length);
        if(nodeIdx < ((*itLines)->length - 1))
          {
            nodeIdx++; // next node
          }
        else
          {
            nodeIdx = 0;
            itLines++; // next line
            lineIdx++;
          }
      }    
    else
    if(itNodes != nodes->end())
      {
        itNodes++;
      }    
  }

  node get()
  {
    node n;
    if(itLines != lines->end())
      {
        n.x = (*itLines)->coord.x + nodeIdx;
        n.y = (*itLines)->coord.y;
        n.z = (*itLines)->coord.z;
      }
    else
      {
        if(itNodes != nodes->end())
          {
            n.x = (*itNodes)->x;
            n.y = (*itNodes)->y;
            n.z = (*itNodes)->z;
          }    
        else
          {
            n.x = 0; n.y = 0; n.z = 0;
          }
      }
    return n;
  }

  void get(node& n)
  {
    if(itLines != lines->end())
      {
        n.x = (*itLines)->coord.x + nodeIdx;
        n.y = (*itLines)->coord.y;
        n.z = (*itLines)->coord.z;
      }
    else
      {
        if(itNodes != nodes->end())
          {
            n.x = (*itNodes)->x;
            n.y = (*itNodes)->y;
            n.z = (*itNodes)->z;
          }    
        else
          {
            n.x = 0; n.y = 0; n.z = 0;
          }
      }
  }

  /*
  void test()
  {
    int nodeSize = 0;
    // count number of pixels in line containers
    for(vector<lineContainer*>::iterator it = lines.begin();
        it != lines.end(); it++)
      nodeSize += (*it)->length;

    printf("nodeSize %d\n", nodeSize);

    int i = 1;
    goToBegin();
    while(!isAtEnd())
      {
        printf("%d ", i);
        i++;
        next();
      }
    printf("\n");
  }
  */

  uint size()
  {
    int nodeSize = 0;
    // count number of pixels in line containers
    for(vector<lineContainer*>::iterator it = lines->begin();
        it != lines->end(); it++)
      nodeSize += (*it)->length;

    return nodeSize + nodes->size();
  }

 private:
  vector<lineContainer*>* lines;
  vector<node*>* nodes;

  // iterators
  vector<lineContainer*>::iterator itLines;
  uint lineIdx;
  vector<node*>::iterator itNodes;
  uint nodeIdx;
};

//------------------------------------------------------------------------------

double node_square_distance(const node& a, const node& b);

//------------------------------------------------------------------------------

class supernode;

/**
 * Supernode aggregates a set of nodes
 */
class supernode
{
 public:
  vector<supernode*> neighbors;
  sidType id;

  // data associated to the supernode.
  supernode_data* data;

  uchar getLabel();

  void setLabel(uchar _label);

  //void setLabel(uchar _label) { if(data) data->label = _label;}

  void setData(uchar _label, int _nr_class, probType* _prob_estimates = 0);

  void addNeighbor(supernode* n)
  {
    neighbors.push_back(n);
  }

  void addLine(lineContainer* l)
  {
    lines.push_back(l);
  }

  void addNode(node* n)
  {
    nodes.push_back(n);
  }

  /**
   * Get center of the supernode with corresponding given id
   * @param center will be initialized by this function
   */
  void getCenter(node& center);

  int getNumberOfNodes() { return nodes.size(); }

  supernode()
  {
    data = 0;
  }

  uint size();

  ~supernode()
    {
      if(data)
        delete data;

      for(vector<node*>::iterator it = nodes.begin();
          it != nodes.end(); it++)
        delete *it;

      for(vector<lineContainer*>::iterator it = lines.begin();
          it != lines.end(); it++)
        delete *it;
    }

  /**
   * Returns a new node iterator
   */
  nodeIterator getIterator() { return nodeIterator(&lines, &nodes); }

 private:
  vector<lineContainer*> lines;
  vector<node*> nodes;

};

#endif // SUPERNODE_H
