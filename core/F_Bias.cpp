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

#include "F_Bias.h"
#include "Config.h"

F_Bias::F_Bias()
{
}

F_Bias::~F_Bias()
{
}

int F_Bias::getSizeFeatureVectorForOneSupernode()
{
  return 1;
}

bool F_Bias::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice, int supernodeId)
{
  x[0].value = 1;
  return true;
}

bool F_Bias::getFeatureVectorForOneSupernode(osvm_node *x, Slice3d* slice3d, int supernodeId)
{
  x[0].value = 1;
  return true;
}

