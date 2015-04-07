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

#include "F_Gaussian.h"
#include "Config.h"

F_Gaussian::F_Gaussian()
{
}

int F_Gaussian::getSizeFeatureVector()
{
  return 2; // mean + variance
}

bool F_Gaussian::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice, int supernodeId)
{
  supernode* s = slice->getSupernode(supernodeId);
  int value;
  node n;

  // Compute mean intensity inside supervoxel
  int nImages = 0;
  int mean = 0;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();

    value = slice->getIntensity(n.x, n.y, n.z);
    // New mean M(n+1) = (n*M(n) + M(n+1))/(n+1)
    mean = (nImages*mean + value)/(nImages+1);
    nImages++;
  }

  // Compute variance
  int temp;
  double variance = 0;
  nImages = 0;
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();

    value = slice->getIntensity(n.x, n.y, n.z);

    // New variance V(n+1) = (n*V(n) + (x(n+1)-M(n))^2)/(n+1)
    temp = value - mean;
    variance = (nImages*variance + temp*temp)/(nImages+1);
    nImages++;
  }

  x[0].value = mean;
  x[1].value = variance;

  return true;
}

bool F_Gaussian::getFeatureVectorForOneSupernode(osvm_node *x, Slice3d* slice3d, int supernodeId)
{
  supernode* s = slice3d->getSupernode(supernodeId);
  int value;
  const ulong size_slice = slice3d->width * slice3d->height;
  node n;

  // Compute mean intensity inside supervoxel
  int nImages = 0;
  int mean = 0;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd())
    {
      ni.get(n);
      ni.next();

      value = slice3d->raw_data[n.z*size_slice + n.y*slice3d->width + n.x];
      // New mean M(n+1) = (n*M(n) + M(n+1))/(n+1)
      mean = (nImages*mean + value)/(nImages+1);
      nImages++;
    }

  // Compute variance
  int temp;
  double variance = 0;
  nImages = 0;
  ni.goToBegin();
  while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();

      value = slice3d->raw_data[n.z*size_slice + n.y*slice3d->width + n.x];

      // New variance V(n+1) = (n*V(n) + (x(n+1)-M(n))^2)/(n+1)
      temp = value - mean;
      variance = (nImages*variance + temp*temp)/(nImages+1);
      nImages++;
    }

  x[0].value = mean;
  x[1].value = variance;

  return true;
}
