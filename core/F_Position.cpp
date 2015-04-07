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

#include "F_Position.h"
#include "Config.h"

F_Position::F_Position(const char* image_name)
{
  if(image_name)
    img = cvLoadImage(image_name);
  else
    img = 0;
}

F_Position::~F_Position()
{
  if(img)
    cvReleaseImage(&img);
}

int F_Position::getSizeFeatureVector()
{
  return 4;
}

bool F_Position::getFeatureVector(osvm_node *n,
                              const int x,
                              const int y)
{
  assert(0);
  n[0].value = (img->width/2)-x;
  n[1].value = (img->height/2)-y;
  return true;
}

bool F_Position::getFeatureVector(osvm_node *x, Slice* slice, int supernodeId)
{
  node center;
  slice->getCenter(supernodeId,center);
  double _x = (center.x - (slice->img_width/2.0))/(double)slice->img_width;
  double _y = center.y/(double)slice->img_height;
  x[0].value = _x;
  x[1].value = _y;
  x[2].value = abs(_x);
  x[3].value = abs(_y);
  return true;
}

bool F_Position::getFeatureVector(osvm_node *x, Slice3d* slice3d, int supernodeId)
{
  assert(0);
  return true;
}

bool F_Position::getFeatureVector(osvm_node *x,
                                  Slice3d* slice3d,
                                  const int gx,
                                  const int gy,
                                  const int gz)
{
  assert(0);
  return true;
}
