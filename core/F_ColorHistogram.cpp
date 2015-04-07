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

#include "F_ColorHistogram.h"
#include "Config.h"

//-------------------------------------------------------------------------TYPES


//---------------------------------------------------------------------FUNCTIONS

F_ColorHistogram::F_ColorHistogram(int _nb_bins,
                                   int _max_pixel_value,
                                   eHistogramType _histoType,
                                   bool _useColorImage,
                                   IplImage* img)
{
  nBinsPerSupernode = _nb_bins;
  nBinsPerSupernode2 = nBinsPerSupernode*nBinsPerSupernode;
  max_pixel_value = _max_pixel_value;
  histoType = _histoType;
  histoType = NO_NEIGHBORS; // temporary hack

  nLocations = 8;
  string config_tmp;
  if(Config::Instance()->getParameter("histogram_nlocations", config_tmp)) {
    nLocations = atoi(config_tmp.c_str());
  }

  // hue and saturation are binned together and value is binned separately
  nBins = nBinsPerSupernode2 + nBinsPerSupernode;

  switch(histoType) {
  case INCLUDE_NEIGHBORS_IN_SEPARATE_BINS:
    {
      offsetNeighbors = nBins;
      nTotalBins = nBins*2;
    }
    break;
  case INCLUDE_NEIGHBORS_PLUS_LOCATION:
    {
      // use nBins*N_LOCATIONS for neighbors
      offsetNeighbors = nBins;
      nTotalBins = nBins + (nBins*nLocations);
    }
    break;
  default:
    {
      offsetNeighbors = 0;
      nTotalBins = nBins;
    }
    break;
  }

  if(img) {
    if(img->nChannels != 3) {
      printf("[F_ColorHistogram] Error : img->nChannels != 3\n");
      exit(-1);
    }
  }

  valToBin = nBinsPerSupernode/(float)max_pixel_value;
}

int F_ColorHistogram::getSizeFeatureVectorForOneSupernode()
{
  return nTotalBins;
}

bool F_ColorHistogram::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice, int supernodeId)
{
  IplImage* _img = slice->colorImg; //FIXME : should probably use a parameter in the config file

  Histogram hist(nTotalBins);

  supernode* s = slice->getSupernode(supernodeId);
  supernode* sn;
  int pValue;
  int idx;
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();

    idx = 0;
    // hue and saturation are binned together
    for(int c=0;c<2;c++) {
      pValue = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+c]);
      idx = idx*nBinsPerSupernode + (pValue*valToBin);
    }
    hist.histData[idx]++;

    // value is binned separately
    pValue = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+2]);
    idx = nBinsPerSupernode2 + (pValue*valToBin);
    hist.histData[idx]++;
  }

  // normalize histogram
  for(int i = 0; i < nBinsPerSupernode2; i++) {
    hist.histData[i] /= (double)s->size();
  }
  for(int i = 0; i < nBinsPerSupernode; i++) {
    hist.histData[nBinsPerSupernode2+i] /= (double)s->size();  
  }

  switch(histoType) {
  case INCLUDE_NEIGHBORS:
    {
      for(vector < supernode* >::iterator itN = s->neighbors.begin();
          itN != s->neighbors.end();itN++) {
        sn = *itN;
        double binValue = 1.0/(double)(2*s->neighbors.size()*sn->size());
        nodeIterator ni = sn->getIterator();
        ni.goToBegin();
        while(!ni.isAtEnd()) {
          ni.get(n);
          ni.next();

          idx = 0;
          // hue and saturation are binned together
          for(int c=0;c<2;c++) {
            pValue = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+c]);
            idx = idx*nBinsPerSupernode + (pValue*valToBin);
          }
          hist.histData[idx]+=binValue;

          // value is binned separately
          pValue = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+2]);
          idx = nBinsPerSupernode2 + (pValue*valToBin);
          hist.histData[idx]+=binValue;
        }
      }
    }
    break;
  case INCLUDE_NEIGHBORS_IN_SEPARATE_BINS:
    {
      for(vector < supernode* >::iterator itN = s->neighbors.begin();
          itN != s->neighbors.end();itN++) {
        sn = *itN;
        double binValue = 1.0/(double)(2*s->neighbors.size()*sn->size());
        nodeIterator ni = sn->getIterator();
        ni.goToBegin();
        while(!ni.isAtEnd()) {
          ni.get(n);
          ni.next();
            
          idx = 0;
          // hue and saturation are binned together
          for(int c=0;c<2;c++) {
            pValue = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+c]);
            idx = idx*nBinsPerSupernode + (pValue*valToBin);
          }
          hist.histData[offsetNeighbors+idx]+=binValue;

          // value is binned separately
          pValue = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+2]);
          idx = nBinsPerSupernode2 + (pValue*valToBin);
          hist.histData[offsetNeighbors+idx]+=binValue;
        }
      }
    }
    break;
  case INCLUDE_NEIGHBORS_PLUS_LOCATION:
    {
      node cs;
      node cn;
      float v[2];
      int angleXY;
      int angleIdx;
      int offsetAngle;
      s->getCenter(cs);

      for(vector < supernode* >::iterator itN = s->neighbors.begin();
          itN != s->neighbors.end();itN++) {
        sn = *itN;

        sn->getCenter(cn);
        // compute vector sn
        v[0] = cs.x - cn.x;
        v[1] = cs.y - cn.y;
        //v[2] = cs.z - cn.z;
                    
        angleXY = atan2(v[1], v[0])*180.0/PI;
        angleXY += 180.0;
        //angleIdx = angleToIdx(angleXY);
        angleIdx = angleXY*nLocations/360.0;
        if(angleIdx >= nLocations) {
          angleIdx = nLocations - 1;
        }
        offsetAngle = offsetNeighbors + (angleIdx*nBins);

        double binValue = 1.0/(double)((nLocations+1)*s->neighbors.size()*sn->size());
        nodeIterator ni = sn->getIterator();
        ni.goToBegin();
        while(!ni.isAtEnd()) {
          ni.get(n);
          ni.next();
            
          idx = 0;
          // hue and saturation are binned together
          for(int c=0;c<2;c++) {
            pValue = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+c]);
            idx = idx*nBinsPerSupernode + (pValue*valToBin);
          }
          hist.histData[offsetAngle+idx]+=binValue;

          // value is binned separately
          pValue = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+2]);
          idx = nBinsPerSupernode2 + (pValue*valToBin);
          hist.histData[offsetAngle+idx]+=binValue;
        }
      }
    }
    break;
  default:
    break;
  }

  for(int i = 0;i < hist.nBins; i++) {
    x[i].value = hist.histData[i];
  }
  
  return true;
}

bool F_ColorHistogram::getFeatureVectorForOneSupernode(osvm_node *x, Slice3d* slice3d, int supernodeId)
{
  printf("[F_ColorHistogram] Error : getFeatureVector(osvm_node *x, Slice3d& slice3d, int supernodeId) not implemented\n");
  exit(-1);

  return false;
}
