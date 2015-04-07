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

#ifndef OSVM_H
#define OSVM_H

#include <stdio.h>
#include <errno.h>
#include "globalsE.h"
#include "Feature.h"
#include "Slice_P.h"
#include "Slice.h"
#include "Slice3d.h"
#include <map>
#include <set>
#include <vector>

#include "svm.h"

#include "oSVM_types.h"

using namespace std;

//---------------------------------------------------------------------- GLOBALS

#define L1_NORM 1
#define L2_NORM 2

//------------------------------------------------------------------------------

class oSVM
{
 public:

  static void compute_range(Slice_P* slice, Feature* feature, map<int, double>& feature_min, map<int, double>& feature_max);

  static void rescale(osvm_node *x, map<int, double>& feature_min, map<int, double>& feature_max);

  static void initSVMNode(osvm_node*& x, int d);

  static void initSVMNodes(osvm_node**& xs, uint nNodes, int d);

  static double norm(osvm_node* x, int norm_type);

  static void print(osvm_node *x, const char* title=0);

  static void printNonZeros(osvm_node *x, const char* title=0);

};

#endif //OSVM_H
