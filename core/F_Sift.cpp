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

#include "F_Sift.h"
#include "highgui.h"

F_Sift::F_Sift(const char* image_name,
               int aoctaves,
               int alevels,
               int aomin)
{
  octaves = aoctaves;
  levels = alevels;
  omin = aomin;
  //scale = omin; // default scale. Use setScale to change it.

  // Extract image data and convert it to floats
  // see std::istream& extractPgm(std::istream& in, PgmBuffer& buffer) in sift.cpp
  IplImage* img = cvLoadImage(image_name,0);
  im_pt = new VL::pixel_t[img->width*img->height];
  uchar pValue;
  int idx = 0;
  for(int y=0;y<img->height;y++)
    for(int x=0;x<img->width;x++)
      {
        pValue = (((uchar*)(img->imageData + y*img->widthStep))[x*img->nChannels]);
        im_pt[idx] = pValue/255.0f; 
        idx++;
      }

  float sigman = .5;
  //sigma0 = 1.6 * powf(2.0f, 1.0f / levels);
  // sigma0 is the value used to smooth the first image at the bottom of the pyramid.
  // sigma0=1.2 is the value given in the SIFT paper
  sigma0 = 1.2f;

  PRINT_MESSAGE("[F_Sift] w=%d h=%d s=%f %f octaves=%d levels=%d omin=%d\n",
                img->width, img->height, 
                sigman, sigma0,
                octaves, levels,
                omin);
  
  for(int o = 0;o<octaves;o++)
    {
      for(int l = 0;l<levels;l++)
        {
          float scale = sigma0*powf(2.0f,o+l/(float)levels);
          scales.push_back(scale);

          PRINT_MESSAGE("[F_Sift] scale=%f\n", scale);
        }
    }

  sift = new VL::Sift(im_pt, img->width, img->height, 
                      sigman, sigma0,
                      octaves, levels,
                      omin, -1, levels+1);

  cvReleaseImage(&img);

  descr_pt = new VL::float_t[F_Sift::desc_size];
}

F_Sift::~F_Sift()
{
  delete[] im_pt;
  delete[] descr_pt;
}

void F_Sift::setScales(vector<float>& _scales)
{
  scales = _scales;
}

void F_Sift::addScale(float _scale)
{
  scales.push_back(_scale);
}

bool F_Sift::getFeatureVectorForOneSupernode(osvm_node *n, Slice* slice, int supernodeId)
{
  VL::float_t angle = 0;

  supernode* s = slice->getSupernode(supernodeId);

  /*
  vector<node*>* v = &(s->nodes);
  int n = v->size();
  int pValue;
  int idx;
  for(vector < node* >::iterator iPt = v->begin();
      iPt != v->end();iPt++)
    {
      pValue = (int)(((uchar*)(slice.img->imageData + (*iPt)->y*slice.img->widthStep))[(*iPt)->x*slice.img->nChannels]);
      idx = (pValue*hist.nBins)/max_pixel_value;
      hist.histData[idx]++;
    }
  */

  node center;
  s->getCenter(center);

  int offset = 0;
  for(vector<float>::iterator itScale = scales.begin();
      itScale != scales.end(); itScale++)
    {
      VL::Sift::Keypoint k = sift->getKeypoint(center.x,center.y,*itScale);

      //printf("computeKeypointDescriptor %d\n", F_Sift::desc_size);
      sift->computeKeypointDescriptor(descr_pt, k, angle);

      for(int i=0;i<desc_size;i++)
        n[offset+i].value = descr_pt[i];

      offset += desc_size;
    }

  return true;
}

bool F_Sift::getFeatureVectorForOneSupernode(osvm_node *n, Slice* slice,
                                             const int x, const int y)
{
  VL::float_t angle = 0;

  int offset = 0;
  for(vector<float>::iterator itScale = scales.begin();
      itScale != scales.end(); itScale++)
    {
      VL::Sift::Keypoint k = sift->getKeypoint(x,y,*itScale);
      sift->computeKeypointDescriptor(descr_pt, k, angle);

      for(int i=0;i<desc_size;i++)
        n[offset+i].value = descr_pt[i];

      offset += desc_size;
    }

  return true;
}

int F_Sift::getSizeFeatureVectorForOneSupernode()
{
  return desc_size*scales.size();
}
