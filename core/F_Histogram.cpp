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

#include "F_Histogram.h"
#include "Config.h"

F_Histogram::F_Histogram(int _nb_bins,
                         int _max_pixel_value,
                         eHistogramType _histoType,
                         bool _useColorImage,
                         IplImage* img)
{
  nBinsPerSupernode = _nb_bins;
  max_pixel_value = _max_pixel_value;
  histoType = _histoType;

  // temporary hack. Remove code related to INCLUDE_NEIGHBORS_IN_SEPARATE_BINS.
  histoType = NO_NEIGHBORS;

  if(histoType == INCLUDE_NEIGHBORS_IN_SEPARATE_BINS) {
    nBins = nBinsPerSupernode*2;
  } else {
    nBins = nBinsPerSupernode;
  }

  string config_tmp;
  useColorImage = false;
  if(Config::Instance()->getParameter("use_color_image", config_tmp)) {
    useColorImage = config_tmp.c_str()[0] == '1';
  }
  else
    useColorImage = _useColorImage;

  if(useColorImage)
    nChannels = 3;
  else
    nChannels = 1;

  normalize_l1_norm = false;
  string param_hist_l1 = Config::Instance()->parameters["hist_normalize_l1_norm"];
  if(param_hist_l1 == "1") {
    PRINT_MESSAGE("[F_Histogram] hist_normalize_l1_norm is set to true\n");
    normalize_l1_norm = true;
  }

  normalize_l2_norm = false;
  string param_hist_l2 = Config::Instance()->parameters["hist_normalize_l2_norm"];
  if(param_hist_l2 == "1") {
    PRINT_MESSAGE("[F_Histogram] hist_normalize_l2_norm is set to true\n");
    normalize_l2_norm = true;
  }

  assert(!(normalize_l1_norm && normalize_l2_norm));
}

int F_Histogram::getSizeFeatureVectorForOneSupernode()
{
  return nBins;
}

bool F_Histogram::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice, int supernodeId)
{
  IplImage* _img;
  if(useColorImage)
    _img = slice->colorImg; //FIXME : should probably use a parameter in the config file
  else
    _img = slice->img;

  if(useColorImage && (_img->nChannels != 3)) {
    printf("[F_Histogram] Error : histogram feature initialized with color information but given image does not have 3 channels\n");
    exit(-1);
  }

  Histogram hist(nBins);

  float valToBin = nBinsPerSupernode/(float)max_pixel_value;
  supernode* s = slice->getSupernode(supernodeId);
  supernode* sn;
  int value;
  int idx;
  node n;
  int offset = 0;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    offset = 0;
    for(int c=0;c<nChannels;c++) {
      value = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+c]);
      idx = value*valToBin;
      hist.histData[offset+idx]++;
      offset += nBinsPerSupernode;
    }
  }

  double binValue = 1.0/(double)s->neighbors.size();
  if(histoType == INCLUDE_NEIGHBORS) {
    for(vector < supernode* >::iterator itN = s->neighbors.begin();
        itN != s->neighbors.end();itN++) {
      sn = *itN;
      nodeIterator ni = sn->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd()) {
        ni.get(n);
        ni.next();
        offset = 0;
        for(int c=0;c<nChannels;c++) {
          value = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+c]);
          idx = value*valToBin;
          hist.histData[offset+idx]+=binValue;
          offset += nBinsPerSupernode;
        }
      }
    }
  }
  else if(histoType == INCLUDE_NEIGHBORS_IN_SEPARATE_BINS) {
    for(vector < supernode* >::iterator itN = s->neighbors.begin();
        itN != s->neighbors.end();itN++) {
      sn = *itN;
      nodeIterator ni = sn->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd()) {
        ni.get(n);
        ni.next();

        offset = nBinsPerSupernode*nChannels;
        for(int c=0;c<nChannels;c++) {
          value = (int)(((uchar*)(_img->imageData + n.y*_img->widthStep))[n.x*_img->nChannels+c]);
          idx = value*valToBin;
          hist.histData[offset+idx]+=binValue;
          offset += nBinsPerSupernode;
        }
      }
    }
  }

  // normalize
  if(normalize_l1_norm) {
    double l1_norm = 0;
    for(int i = 0;i < nBinsPerSupernode; i++) {
      l1_norm += hist.histData[i];
    }
    //if(l1_norm != 0) {
    if(fabs(l1_norm) > 1e-20) {
      for(int i = 0;i < nBinsPerSupernode; i++) {
        hist.histData[i] /= l1_norm;
      }
    }
    if(histoType == INCLUDE_NEIGHBORS_IN_SEPARATE_BINS) {
      l1_norm = 0;
      for(int i = 0;i < nBinsPerSupernode; i++) {
        l1_norm += hist.histData[nBinsPerSupernode+i];
      }
      if(fabs(l1_norm) > 1e-20) {
        //if(l1_norm != 0) {
        for(int i = 0;i < nBinsPerSupernode; i++) {
          hist.histData[nBinsPerSupernode+i] /= l1_norm;
        }
      }
    }
  }

  if(normalize_l2_norm) {
    double l2_norm = 0;
    for(int i = 0;i < nBinsPerSupernode; i++) {
      l2_norm += hist.histData[i]*hist.histData[i];
    }
    if(fabs(l2_norm) > 1e-20) {
      for(int i = 0;i < nBinsPerSupernode; i++) {
        hist.histData[i] /= l2_norm;
      }
    }
    if(histoType == INCLUDE_NEIGHBORS_IN_SEPARATE_BINS) {
      l2_norm = 0;
      for(int i = 0;i < nBinsPerSupernode; i++) {
        l2_norm += hist.histData[nBinsPerSupernode+i] * hist.histData[nBinsPerSupernode+i];
      }
      if(fabs(l2_norm) > 1e-20) {
        for(int i = 0;i < nBinsPerSupernode; i++) {
          hist.histData[nBinsPerSupernode+i] /= l2_norm;
        }
      }
    }
  }

  for(int i = 0;i < hist.nBins; i++) {
    x[i].value = hist.histData[i];
  }

  return true;
}

bool F_Histogram::getFeatureVectorForOneSupernode(osvm_node *x, Slice3d* slice3d, int supernodeId)
{
  if(useColorImage && (slice3d->nChannels != 3)) {
    printf("[F_Histogram] Error : histogram feature initialized with color information but given image does not have 3 channels\n");
    exit(-1);
  }

  Histogram hist(nBins*nChannels);

  double valToIdx;
  if(histoType == INCLUDE_NEIGHBORS_IN_SEPARATE_BINS) {
    valToIdx = (hist.nBins/2.0)/(double)max_pixel_value;
  } else {
    valToIdx = hist.nBins/(double)max_pixel_value;
  }

  uchar* raw_data = slice3d->getRawData();
  supernode* s = slice3d->getSupernode(supernodeId);
  supernode* sn;
  int value;
  int idx;
  const ulong size_slice = slice3d->width * slice3d->height;
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    for(int c=0;c<nChannels;c++) {
      value = raw_data[n.z*size_slice + n.y*slice3d->width + n.x + c];
      idx = value*valToIdx;
      hist.histData[idx]++;
    }
  }

  double binValue = 1.0/(double)s->neighbors.size();
  if(histoType == INCLUDE_NEIGHBORS) {
    for(vector < supernode* >::iterator itN = s->neighbors.begin();
        itN != s->neighbors.end();itN++) {
      sn = *itN;
      nodeIterator ni = sn->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd()) {
        ni.get(n);
        ni.next();
        
        for(int c=0;c<nChannels;c++) {
          value = raw_data[n.z*size_slice + n.y*slice3d->width + n.x + c];
          idx = value*valToIdx;
          hist.histData[idx]+=binValue;
        }
      }
    }
  } else{
    if(histoType == INCLUDE_NEIGHBORS_IN_SEPARATE_BINS) {
      for(vector < supernode* >::iterator itN = s->neighbors.begin();
          itN != s->neighbors.end();itN++) {
        sn = *itN;
        nodeIterator ni = sn->getIterator();
        ni.goToBegin();
        while(!ni.isAtEnd()) {
          ni.get(n);
          ni.next();
          
          for(int c=0;c<nChannels;c++) {
            value = raw_data[n.z*size_slice + n.y*slice3d->width + n.x + c];
            idx = value*valToIdx;
            hist.histData[nBinsPerSupernode+idx]+=binValue;
          }
        }
      }
    }
  }

  // normalize
  if(normalize_l1_norm) {
    double l1_norm = 0;
    for(int i = 0;i < nBinsPerSupernode; i++) {
      l1_norm += hist.histData[i];
    }
    for(int i = 0;i < nBinsPerSupernode; i++) {
       hist.histData[i] /= l1_norm;
    }
    if(histoType == INCLUDE_NEIGHBORS_IN_SEPARATE_BINS) {
      l1_norm = 0;
      for(int i = 0;i < nBinsPerSupernode; i++) {
        l1_norm += hist.histData[nBinsPerSupernode+i];
      }
      if(fabs(l1_norm) > 1e-20) {
        for(int i = 0;i < nBinsPerSupernode; i++) {
          hist.histData[nBinsPerSupernode+i] /= l1_norm;
        }
      }
    }
  }

  if(normalize_l2_norm) {
    double l2_norm = 0;
    for(int i = 0;i < nBinsPerSupernode; i++) {
      l2_norm += hist.histData[i]*hist.histData[i];
    }
    if(fabs(l2_norm) > 1e-20) {
      for(int i = 0;i < nBinsPerSupernode; i++) {
        hist.histData[i] /= l2_norm;
      }
    }
    if(histoType == INCLUDE_NEIGHBORS_IN_SEPARATE_BINS) {
      l2_norm = 0;
      for(int i = 0;i < nBinsPerSupernode; i++) {
        l2_norm += hist.histData[nBinsPerSupernode+i] * hist.histData[nBinsPerSupernode+i];
      }
      if(fabs(l2_norm) > 1e-20) {
        for(int i = 0;i < nBinsPerSupernode; i++) {
          hist.histData[nBinsPerSupernode+i] /= l2_norm;
        }
      }
    }
  }

  for(int i = 0;i < hist.nBins; i++) {
    x[i].value = hist.histData[i];
  }

  return true;
}

