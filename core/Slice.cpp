////////////////////////////////////////////////////////////////////////
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

#include "Slice.h"
#include "utils.h"

#include "LKM.h"

// standard libraries
#include <limits.h>
#include <sstream>
#include <time.h>

//------------------------------------------------------------------------------

#define USE_UNDIRECTED_GRAPH 0

//------------------------------------------------------------------------------

void Slice::init(int width, int height)
{
  PRINT_MESSAGE("[Slice] Initializing slice. width %d height %d\n",width,height);
  img_width = width;
  img_height = height;
  min_sid = INT_MAX; //numeric_limits<int>::max();

  // Populating mSupernodes
  map<sidType, supernode* >::iterator iLabel;
  int iBuffer = 0;
  sidType sid;  // supernode id
  for(int y=0;y<img_height;y++) {
    for(int x=0;x<img_width;x++) {
      node* p = new node;
      p->x = x;
      p->y = y;
      p->z = 0;
      
      sid = pixelLabels[iBuffer];
      iLabel = mSupernodes.find(sid);
      if(iLabel != mSupernodes.end()) {
        iLabel->second->addNode(p);
      } else {
        // Create a new supernode
        supernode* s = new supernode;
        s->id = sid;
        s->addNode(p);
        mSupernodes[sid] = s;

        if(s->id < min_sid)
          min_sid = s->id;
      }

      iBuffer++;
    }
  }

  if(min_sid > 0) {
    printf("[Slice] WARNING : min_sid equals %d. Should be 0 ?\n", min_sid);
  }

  // initialize random seed
  srand(time(0));

  // edges are set by loadNeighborhoodMap
  nbEdges = 0;
  avgIntensity = -1;
  supernodeLabelsLoaded = false;
  neighborhoodMapLoaded = false;

  minPercentToAssignLabel = MIN_PERCENT_TO_ASSIGN_LABEL;
  includeOtherLabel = minPercentToAssignLabel != 0;

  // load a color image as well
  generateColorImage();
}

Slice::Slice(const char* fn_label, int awidth, int aheight)
{
  img = 0;
  colorImg = 0;

  // Loading file containing a superpixel label for each pixel
  int bufferSize = awidth*aheight;
  pixelLabels = new sidType[bufferSize];
  ifstream ifslb(fn_label, ios::in | ios::binary);
  if(ifslb.fail()) {
    printf("[Slice] Error while loading %s\n",fn_label);
    return;
  }
  ifslb.read((char*)pixelLabels, bufferSize*sizeof(sidType));
  ifslb.close();

  init(awidth,aheight);
}

Slice::Slice(const char* a_image_name,
             const char* fn_label,
             int superpixelStepSize,
             float M)
{
  initSuperpixels(a_image_name, fn_label, superpixelStepSize, M, true);
}

Slice::Slice(const char* a_image_name, const char* fn_label,
             int superpixelStepSize, float M, bool bGenerateNeighborhoodMap)
{
  initSuperpixels(a_image_name, fn_label, superpixelStepSize, M,
                  bGenerateNeighborhoodMap);
}

void Slice::initSuperpixels(const char* a_image_name, const char* fn_label,
                            int superpixelStepSize, float M,
                            bool bGenerateNeighborhoodMap)
{
  colorImg = 0;
  supernode_step = superpixelStepSize;
  cubeness = M;

  bool superpixelLoaded = false;
  if(fn_label != 0 && fileExists(fn_label)) {
    image_name = (string)a_image_name;
    img = cvLoadImage(a_image_name,CV_LOAD_IMAGE_COLOR);
    eraseImage = true;
    if(!img) {
      printf("[Slice] Error : input image %s was not found\n", a_image_name);
      return;
    }

    // Loading file containing a superpixel label for each pixel
    int bufferSize = img->width*img->height;
    pixelLabels = new sidType[bufferSize];

    string labelFilename;
    if(isDirectory(fn_label)) {
      labelFilename = fn_label;
      labelFilename += "/";
      labelFilename += getNameFromPathWithoutExtension(a_image_name);
      labelFilename += ".dat";
    } else {
      labelFilename = fn_label;
    }

    PRINT_MESSAGE("[Slice] Loading %s\n", labelFilename.c_str());
    ifstream ifslb(labelFilename.c_str(), ios::in | ios::binary);
    if(ifslb.fail()) {
      printf("[Slice] Error while loading %s\n", labelFilename.c_str());
    } else {
      string imageNameWithoutPath = getNameFromPath(string(a_image_name));
      ifslb.read((char*)pixelLabels, bufferSize*sizeof(sidType));
      ifslb.close();

      init(img->width,img->height);

      if(bGenerateNeighborhoodMap) {
        generateNeighborhoodMap(pixelLabels, img->width, img->height);
        neighborhoodMapLoaded = true;
      }
      superpixelLoaded = true;
    }
  }
  
  if (!superpixelLoaded) {
      generateSuperpixels(a_image_name, superpixelStepSize, M);

      if(bGenerateNeighborhoodMap) {
        generateNeighborhoodMap(pixelLabels, img_width, img_height);
        neighborhoodMapLoaded = true;
      }
    }
}

Slice::Slice(IplImage* _img, const char* fn_label, int superpixelStepSize, float M)
{
  colorImg = 0;
  image_name = "";
  img = _img;
  eraseImage = false;
  supernode_step = superpixelStepSize;
  cubeness = M;

  if(fn_label != 0 && fileExists(fn_label)) {
    // Loading file containing a superpixel label for each pixel
    int bufferSize = img->width*img->height;
    pixelLabels = new sidType[bufferSize];
    ifstream ifslb(fn_label, ios::in | ios::binary);
    if(ifslb.fail()) {
      printf("[Slice] Error while loading %s\n",fn_label);
      return;
    }
    ifslb.read((char*)pixelLabels, bufferSize*sizeof(sidType));
    ifslb.close();
    init(img->width,img->height);
  } else {
    generateSuperpixels(superpixelStepSize, M);
  }

  generateNeighborhoodMap(pixelLabels, img_width, img_height);
  neighborhoodMapLoaded = true;
}

Slice::Slice(const char* a_img_name, int superpixelStepSize, float M)
{
  colorImg = 0;
  eraseImage = false;
  supernode_step = superpixelStepSize;
  cubeness = M;
  generateSuperpixels(a_img_name, superpixelStepSize, M);
  generateNeighborhoodMap(pixelLabels, img_width, img_height);
  neighborhoodMapLoaded = true;
}

void Slice::generateSuperpixels(const char* a_img_name, int superpixelStepSize, float M)
{
  image_name = (string)a_img_name;
  img = cvLoadImage(a_img_name,CV_LOAD_IMAGE_COLOR);
  eraseImage = true;
  if(!img) {
    printf("[Slice] Error : input image %s was not found\n", a_img_name);
    return;
  }
  generateSuperpixels(superpixelStepSize, M);
}

void Slice::generateSuperpixels(int superpixelStepSize, float M)
{
  img_width = img->width;
  img_height = img->height;
  supernode_step = superpixelStepSize;
  cubeness = M;

  switch(cubeness) {
  case SUPERPIXEL_IMAGE:
    {
      uint img_size = img_width*img_height;
      pixelLabels = new sidType[img_size];
      for(int i =0; i < img_size; ++i) {
        pixelLabels[i] = 0;
      }
      init(img->width,img->height);
    }    
    break;
  case SUPERPIXEL_CUBE:
    {
      // create cubes
      pixelLabels = new sidType[img_width*img_height];
      int sid = 0;
      for(int x = 0; x < img_width; x += superpixelStepSize) {
        for(int y = 0; y< img_height; y += superpixelStepSize) {
          for(int sx=x;sx<min(img_width,(sizeSliceType)x+superpixelStepSize);sx++) {
            for(int sy=y;sy<min(img_height,(sizeSliceType)y+superpixelStepSize);sy++) {
              pixelLabels[sy*img_width+sx] = sid;
            }
          }
          sid++;
        }
      }
      init(img->width,img->height);

    }
    break;
  default:
    {
      int n = img->height*img->width;
      uint* ubuff = new uint[n];
      uint pValue = 0;
      char c;
      uint r,g,b;
      int idx = 0;
      for(int j=0;j<img->height;j++)
        for(int i=0;i<img->width;i++)
          {
            if(img->nChannels == 3)
              {
                // image is assumed to have data in BGR order
                b = ((uchar*)(img->imageData + img->widthStep*(j)))[(i)*img->nChannels];
                g = ((uchar*)(img->imageData + img->widthStep*(j)))[(i)*img->nChannels+1];
                r = ((uchar*)(img->imageData + img->widthStep*(j)))[(i)*img->nChannels+2];
                pValue = b | (g << 8) | (r << 16);
              }
            else if(img->nChannels == 1)
              {
                c = ((uchar*)(img->imageData + img->widthStep*(j)))[(i)*img->nChannels];
                pValue = c | (c << 8) | (c << 16);
              }
            else
              {
                printf("[Slice] Unknown number of channels %d\n", img->nChannels);
                assert(0);
              }          
            ubuff[idx] = pValue;
            idx++;
          }

      LKM* lkm = new LKM;
      int nSupernodes;
      PRINT_MESSAGE("[Slice] Generating SLIC superpixels. superpixelStepSize=%d, M=%f\n", superpixelStepSize, M);
      lkm->DoSuperpixelSegmentation(ubuff,
                                    (const int)img_width, (const int)img_height,
                                    pixelLabels, nSupernodes,
                                    superpixelStepSize, M);

      delete lkm;
      init(img->width,img->height);
      delete[] ubuff;
    }
    break;
  }
}

Slice::~Slice()
{
  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
    delete it->second;
  }

  delete[] pixelLabels;

  if(img != 0) {
    if(eraseImage == false) {
      PRINT_MESSAGE("[Slice] Do not erase image\n");
    } else {
      cvReleaseImage(&img);
    }
  }

  if(colorImg != 0) {
    cvReleaseImage(&colorImg);
  }
}

void Slice::exportSuperpixels(const char* filename)
{
  ofstream ofs(filename, ios::binary);
  ofs.write((char*)(pixelLabels), sizeof(sidType)*img->width*img->height);
  ofs.close();
}

bool Slice::loadNeighborhoodMap(const char* fn_neighbors)
{
  if(!neighborhoodMapLoaded)
    {
      PRINT_MESSAGE("[Slice] Loading neighborhood map %s\n", fn_neighbors);
      const int MAX_LENGTH = 2000;
      char line[MAX_LENGTH];
      int key; // supernode key
      int n; // neighbors key

      ifstream ifs(fn_neighbors);
      if(!ifs) {
        printf("[Slice] Error while loading %s\n",fn_neighbors);
        printf("[Slice] Generating neighborhood map\n");
        generateNeighborhoodMap(pixelLabels, img_width, img_height);
      } else {
          nbEdges = 0;
          while(ifs.getline(line, MAX_LENGTH)) {
            istringstream iss(line);
            iss >> key;
            supernode* s = getSupernode(key);

            // load neighbors
            while(!(iss >> n).fail()) {
              supernode* sn = getSupernode(n);
              s->neighbors.push_back(sn);
              nbEdges++;
            }
          }
          ifs.close();
        }

      nbEdges /= 2; // number of undirected edges
    }

  return true;
}

bool Slice::generateNeighborhoodMap()
{
  if(!neighborhoodMapLoaded) {
    generateNeighborhoodMap(pixelLabels, img->width, img->height);
    neighborhoodMapLoaded = true;
  }
  return neighborhoodMapLoaded;
}

bool Slice::generateNeighborhoodMap(sidType* klabels,
                                    const int width,
                                    const int height)
{

#if USE_LONG_RANGE_EDGES
  addLongRangeEdges();
  return true;
#endif

  int nSupernodes = mSupernodes.size();
  PRINT_MESSAGE("[Slice] Generating neighborhood map for %d labels.\n", nSupernodes);

  const int nh_size = 1; // neighborhood size
  sidType sid;
  sidType nsid;
  vector<supernode*>* ptrNeighbors;
  supernode* s;
  supernode* sn;
  bool existingSupernode;
  nbEdges = 0;
  for(int x = nh_size; x < width - nh_size; x++) {
    for(int y = nh_size; y < height - nh_size; y++) {
      sid = klabels[y*width+x];
      for(int nx = x-nh_size; nx <= x+nh_size; nx++) {
        for(int ny = y-nh_size; ny <= y+nh_size; ny++) {
          nsid = klabels[ny*width+nx];
          if(sid > nsid) {
            s = mSupernodes[sid];
            sn = mSupernodes[nsid];
            if(sn == 0) {
              printf("[Slice3d] Error : supernode %d is null (coordinate=(%d,%d))\n",nsid,nx,ny);
              exit(-1);
            }
            ptrNeighbors = &(s->neighbors);
            existingSupernode = false;
            // Searching in a hash map would be faster but it would also consume
            // more memory ...
            for(vector<supernode*>::iterator itN = ptrNeighbors->begin();
                itN != ptrNeighbors->end(); itN++) {
              if((*itN)->id == nsid) {
                existingSupernode = true;
                break;
              }
            }
            if(!existingSupernode) {
              s->neighbors.push_back(sn);
              sn->neighbors.push_back(s);
              nbEdges++;
            }
          }
        }
      }
    }
  }

  // compute average degree
  double avgDegree = 0;
  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
    s = it->second;
    avgDegree += s->neighbors.size();
  }
  avgDegree /= mSupernodes.size();
  PRINT_MESSAGE("[Slice] %ld supernodes avgDegree %g\n", mSupernodes.size(), avgDegree);

  return true;
}

bool Slice::computeStats(int supernodeId, double& pMean, double& pVar)
{
  supernode* s = getSupernode(supernodeId);
  int pValue;
  pMean = 0;
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    pValue = (int)(((uchar*)(img->imageData + n.y*img->widthStep))[n.x*img->nChannels]);
    pMean += pValue;
  }

  pMean /= s->size();

  // Compute variance
  pVar = 0;
  ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();

    pValue = (int)(((uchar*)(img->imageData + n.y*img->widthStep))[n.x*img->nChannels]) - pMean;
    pVar += pValue*pValue;
  }
  pVar /= s->size();

  return true;
}

bool Slice::getCenter(int supernodeId, node& center)
{
  supernode* s = getSupernode(supernodeId);
  s->getCenter(center);
  return true;
}

ulong Slice::getNbEdges()
{
  return nbEdges;
}

ulong Slice::getNbSupernodes()
{
  return mSupernodes.size();
}

ulong Slice::getNbNodes()
{
  return img_height*img_width;
}

int Slice::getIntensity(int x, int y, int z)
{
  // returns intensity of the first channel
  return (int)cvGet2D(img, y, x).val[0];
}

float Slice::getAvgIntensity(sidType supernodeId)
{
  supernode* s = getSupernode(supernodeId);
  float intensity = 0;
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    for(int c = 0; c < img->nChannels; c++)
      intensity += (float)cvGet2D(img, n.y, n.x).val[c];
  }
  intensity /= s->size()*img->nChannels;
  return intensity;
}

float Slice::getAvgIntensity(int supernodeId, int& r, int &g, int &b)
{
  supernode* s = getSupernode(supernodeId);

  IplImage* _img = img;
  // check if color image was loaded
  if(colorImg) {
    _img = colorImg;
  }

  ulong lr = 0;
  ulong lg = 0;
  ulong lb = 0;

  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    lr += cvGet2D(_img, n.y, n.x).val[2];
    lg += cvGet2D(_img, n.y, n.x).val[1];
    lb += cvGet2D(_img, n.y, n.x).val[0];
  }
  ulong t = s->size();
  r = lr/t;
  g = lg/t;
  b = lb/t;
  return 0;
}

labelType Slice::computeSupernodeLabel(supernode* s, IplImage* imAnnotation)
{
  labelType label;
  int countObject = 0;
  int countBackground = 0;
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    if(((uchar*)(imAnnotation->imageData + n.y*imAnnotation->widthStep))[n.x*imAnnotation->nChannels] == BACKGROUND_MASKVALUE)
      countBackground++;
    else
      countObject++;
  }

  if(includeOtherLabel) {
    int total = ni.size();
    if(countObject>minPercentToAssignLabel*total) {
        label = FOREGROUND;
    } else {
      if(countBackground>minPercentToAssignLabel*total)
        label = BACKGROUND;
      else {
        label = OTHER_LABEL;
      }
    }
  } else {
    if(countObject>countBackground) {
      label = FOREGROUND;
    } else {
      label = BACKGROUND;
    }
  }
  
  return label;
}

void Slice::generateSupernodeLabelFromTextFile(const char* fn_annotation,
                                               int _nLabels)
{
  ifstream ifs(fn_annotation);
  if(!ifs) {
    printf("[Slice] Error while loading text file %s\n",fn_annotation);
    return;
  }

  const int MAX_LENGTH = 100;
  char line[MAX_LENGTH];
  int label;
  supernode* s;
  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
    s = it->second;
    ifs.getline(line, MAX_LENGTH);
    label = (int)atof(line);
    s->setData(label,_nLabels);
  }
  ifs.close();
}

void Slice::generateSupernodeLabels(const char* fn_annotation,
                                    bool includeBoundaryLabels,
                                    bool useColorImages)
{
  if(supernodeLabelsLoaded) {      
    printf("[Slice] Warning : Supernode labels have already been loaded\n");
    return;
  }

  if(useColorImages) {
    // use dummy classIdxToLabel map
    map<ulong, labelType> classIdxToLabel;
    generateSupernodeLabelsFromMultiClassMaskImage(fn_annotation,
                                                   classIdxToLabel);
  } else {
    generateSupernodeLabelFromMaskImage(fn_annotation,
                                        includeBoundaryLabels);
  }
  supernodeLabelsLoaded = true;
}

void Slice::generateSupernodeLabelFromMaskImage(const char* fn_annotation,
                                                bool includeBoundaryLabels)
{
  // Load annotation image
  IplImage* imAnnotation = cvLoadImage(fn_annotation);
  if(!imAnnotation) {
    printf("[Slice] Error while loading %s\n", fn_annotation);
    return;
  }

  supernode* s;
  nLabels = 2;
  if(includeBoundaryLabels) {
    nLabels = 3;
  }

  // first pass to compute labels
  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
    s = it->second;
    s->setData(computeSupernodeLabel(s, imAnnotation),nLabels);
  }

  if(includeBoundaryLabels) {
    // second pass to change foreground labels to boundary labels if they are
    // touching background labels
    for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
        it != mSupernodes.end(); it++) {
      s = it->second;
      if(s->getLabel() == FOREGROUND) {
        for(vector < supernode* >::iterator itN = s->neighbors.begin();
            itN != s->neighbors.end();itN++) {
          supernode* ns = *itN;
          if(s->id == ns->id) {
            printf("[Slice] Error : supernode %d has a neighbor with the same sid\n", s->id);
            continue;
          }

          if(ns->getLabel() == BACKGROUND) {
            // BACKGROUND label was found among the neighbors
            // switch label to BOUNDARY
            s->setLabel(BOUNDARY);
            break;
          }
        }
      }
    }
  }

  cvReleaseImage(&imAnnotation);
}

ulong Slice::computeSupernodeLabelFromMulticlassImage(sidType sid,
                 IplImage* mask, map<ulong, labelType>& classIdxToLabel)
{
  ulong sLabel = 0; // returned label
  ulong classIdx;
  map<ulong, ulong> lPixelsPerClass;
  uchar pixelValue;
  supernode* s = mSupernodes[sid];
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();
    classIdx = 0;
    for(int c = 0; c < mask->nChannels; c++) {
      pixelValue = ((uchar*)(mask->imageData + mask->widthStep*n.y))[n.x*mask->nChannels+c];
      classIdx += pow(255.0,c)*pixelValue;
    }
    
    if(lPixelsPerClass.count(classIdx) == 0) {
      lPixelsPerClass[classIdx] = 0;
    } else {
      lPixelsPerClass[classIdx]++;
    }
  }

  // pick class with maximum number of pixels
  ulong maxPixels = 0;
  bool init = false;
  for(map<ulong, ulong>::iterator itClass = lPixelsPerClass.begin();
      itClass != lPixelsPerClass.end(); itClass++) {
    if(!init || itClass->second > maxPixels) {
      init = true;
      maxPixels = itClass->second;

      if (classIdxToLabel.count(itClass->first) == 0) {
        sLabel = itClass->first;
      } else {
        sLabel = classIdxToLabel[itClass->first];
      }
    }
  }

  return sLabel;
}

void Slice::generateSupernodeLabelsFromMultiClassMaskImage(const char* fn_annotation,
                                                           map<ulong, labelType>& classIdxToLabel)
{
  // Load annotation image
  IplImage* imAnnotation = cvLoadImage(fn_annotation);
  if(!imAnnotation) {
    printf("[Slice] Error while loading %s\n",fn_annotation);
    return;
  }

  // compute labels
  supernode* s;
  labelType label;
  int nLabels = classIdxToLabel.size();
  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
      s = it->second;
      label = computeSupernodeLabelFromMulticlassImage(s->id, imAnnotation, classIdxToLabel);
      s->setData(label, nLabels);
    }

  cvReleaseImage(&imAnnotation);
}


int Slice::search(int x, int y)
{
  float dx,dy,d;
  float minDist = FLT_MAX;
  int sid = -1;
  node n;
  supernode* s;

  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
    s = it->second;
    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();

      dx = n.x - x;
      dy = n.y - y;
      d = dx*dx + dy*dy;
      if(d<minDist) {
        minDist = d;
        sid = it->first;
      }
    }
  }
  return sid;
}

void Slice::generateColorImage()
{
  if(colorImg == 0 && img->nChannels == 3) {
    colorImg = cvCreateImage(cvGetSize(img),IPL_DEPTH_8U,3);

    // OpenCv resizes each channel to be between 0 and 255
    //cvCvtColor(img,colorImg,CV_RGB2Lab);
    cvCvtColor(img,colorImg,CV_BGR2HSV);
  }
}

void Slice::exportTextLabels(const char* filename)
{
  ofstream ofs(filename);
  for(map<sidType, supernode* >::iterator it = this->mSupernodes.begin();
      it != this->mSupernodes.end(); it++) {
    ofs << (int)it->second->getLabel() << endl;
  }
  ofs.close();
}

void Slice::exportSupernodeLabels(const char* filename, int nClasses,
			   labelType* labels,
			   int nLabels,
			   const map<labelType, ulong>* labelToClassIdx)
{
  createColoredAnnotationImage(filename, nClasses, labels, nLabels,
			       labelToClassIdx);
}

int Slice::createColoredAnnotationImage(const char* outputFilename,
                                        int nstates,
                                        labelType* labels,
                                        int nLabels,
                                        const map<labelType, ulong>* labelToClassIdx)
{
  IplImage* img = getColoredAnnotationImage(nstates,
                                            labels,
                                            nLabels,
                                            labelToClassIdx);
  cvSaveImage(outputFilename,img);
  cvReleaseImage(&img);
  return 0;
}

IplImage* Slice::getColoredAnnotationImage(int nstates,
                                           labelType* labels,
                                           int nLabels,
                                           const map<labelType, ulong>* labelToClassIdx)
{
  IplImage* img = cvCreateImage(cvSize(this->img_width, this->img_height),
                                IPL_DEPTH_8U,3);

  supernode* s;
  uchar *pData;
  node n;
  uchar r,g,b;
  ulong classIdx = 0;
  int label;
  int *counter = new int[nstates];
  for(int i = 0; i < nstates; i++) {
    counter[i]=0;
  }

  for(map<sidType, supernode* >::iterator it = this->mSupernodes.begin();
      it != this->mSupernodes.end(); it++) {
    label = labels[it->first];
    counter[label]++;

    if(labelToClassIdx == 0) {
      r = g = b = label;
    } else {
      classIdx = labelToClassIdx->find(label)->second;
      classIdxToRGB(classIdx,r,g,b);
    }

    s = it->second;
    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd())
      {
        ni.get(n);
        ni.next();
        pData = (((uchar*)(img->imageData + n.y*img->widthStep)+n.x*img->nChannels));
        pData[0] = b;
        pData[1] = g;
        pData[2] = r;
      }
  }

  /*
  PRINT_MESSAGE("[Slice] ");
  for(int i = 0; i < nstates; i++)
    PRINT_MESSAGE("local[%d]=%d ",i,counter[i]);
  PRINT_MESSAGE("\n");
  */

  delete [] counter;
  return img;
}

void Slice::exportOverlay(const char* filename)
{
  if(!this->img) {
    printf("[Slice] Error : no image provided\n");
    exit(-1);
  }
  IplImage* _img = cvCloneImage(this->img);

  uchar *pData;
  int label;
  node n;
  supernode* s;
  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
    s = it->second;
    label = s->getLabel();
    if(label == T_FOREGROUND) {
      nodeIterator ni = s->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd()) {
        ni.get(n);
        ni.next();
        
        pData = (((uchar*)(_img->imageData + n.y*_img->widthStep)+n.x*_img->nChannels));
        pData[0] = 255;
        //pData[1] = g;
        //pData[2] = r;
      }
    }
  }

  cvSaveImage(filename, _img);
  cvReleaseImage(&_img);
}

void Slice::exportOverlay(const char* filename, labelType* labels)
{
  if(!this->img) {
    printf("[Slice] Error : no image provided\n");
    exit(-1);
  }
  IplImage* _img = cvCloneImage(this->img);

  uchar *pData;
  int label;
  node n;
  supernode* s;
  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
    label = labels[it->first];
    if(label == T_FOREGROUND) {
      s = it->second;
      nodeIterator ni = s->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd()) {
        ni.get(n);
        ni.next();
        
        pData = (((uchar*)(_img->imageData + n.y*_img->widthStep)+n.x*_img->nChannels));
        pData[0] = 255;
        //pData[1] = g;
        //pData[2] = r;
      }
    }
  }

  cvSaveImage(filename, _img);
  cvReleaseImage(&_img);
}

void Slice::exportProbabilities(const char* filename, int nClasses,
                                float* pbs)
{
  IplImage* _img = cvCreateImage(cvGetSize(img),IPL_DEPTH_8U,1);

  uchar *pData;
  float pb;
  node n;
  supernode* s;
  for(map<sidType, supernode* >::iterator it = mSupernodes.begin();
      it != mSupernodes.end(); it++) {
    pb = pbs[it->first];

    s = it->second;
    nodeIterator ni = s->getIterator();
    ni.goToBegin();
    while(!ni.isAtEnd()) {
      ni.get(n);
      ni.next();
        
      pData = (((uchar*)(_img->imageData + n.y*_img->widthStep)+n.x*_img->nChannels));
      pData[0] = pb*255;
    }
  }

  cvSaveImage(filename, _img);
  cvReleaseImage(&_img);  
}
