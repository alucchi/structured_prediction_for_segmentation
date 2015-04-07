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

#ifndef SLICE_H
#define SLICE_H

// standard libraries
#include <cmath> //std::log...
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cv.h>
#include <highgui.h>

// SliceMe
#include "colormap.h"
#include "globalsE.h"
#include "Slice_P.h"
#include "Supernode.h"

using namespace std;

//------------------------------------------------------------------------------

/**
 * Slice is a class used to stored a Slice made of "supernodes"
 * Each supernode is composed of nodes (pixels).
 */
class Slice : public Slice_P
{
 public:
  IplImage* img; // source image
  int img_width;
  int img_height;
  sidType* pixelLabels;
  string image_name;

  /**
   * Color image : some features need take an image in the LAB or HSV space as input
   */
  IplImage* colorImg;

  // map containing a list of supernodes indexed by their supernode id
  map<sidType, supernode* > mSupernodes;

  Slice() {}

  Slice(const char* a_image_name, const char* fn_label,
        int superpixelStepSize, float M, bool bGenerateNeighborhoodMap);

  Slice(const char* a_img_name,
        const char* fn_label,
        int superpixelStepSize = SUPERPIXEL_DEFAULT_STEP_SIZE,
        float M = SUPERPIXEL_DEFAULT_M);

  /**
   * Load superpixel labels from a given file
   * @param img input image
   * @param fn_label is a file containing the labels for each pixel
   */
  Slice(IplImage* _img,
        const char* fn_label = 0,
        int superpixelStepSize = SUPERPIXEL_DEFAULT_STEP_SIZE,
        float M = SUPERPIXEL_DEFAULT_M);

  /**
   * Load superpixel labels from a given file
   * @param a_width is the width of the image
   * @param a_height is the height of the image
   * @param fn_label is a file containing the labels for each pixel
   */
  Slice(const char* fn_label, int awidth, int aheight);

  /**
   * Generate superpixel labels for the given file
   * @param img_name name of the image
   */
  Slice(const char* img_name,
        int superpixelStepSize = SUPERPIXEL_DEFAULT_STEP_SIZE,
        float M = SUPERPIXEL_DEFAULT_M);

  /**
   * Destructor frees the allocated resources.
   */
  ~Slice();

#if ADD_LONG_EDGES

  /**
   * Add long range edges
   */
  void addLongRangeEdges();
#endif

  bool computeStats(int supernodeId,
                    double& pMean, double& pVar);

  int createColoredAnnotationImage(const char* outputFilename,
                                   int nstates,
                                   labelType* labels,
                                   int nLabels,
                                   const map<labelType, ulong>* labelToClassIdx);

  void exportOverlay(const char* filename);

  void exportOverlay(const char* filename, labelType* labels);

  void exportProbabilities(const char* filename, int nClasses,
                           float* pbs);

  void exportSupernodeLabels(const char* filename, int nClasses,
                             labelType* labels,
                             int nLabels,
			     const map<labelType, ulong>* labelToClassIdx);

  void exportSuperpixels(const char* filename);

  /**
   * Create a text file with the labels of the superpixels.
   * Print one label per line.
   */
  void exportTextLabels(const char* filename);

  void generateColorImage();

  bool generateNeighborhoodMap();

  /**
   * Generate neighborhood map
   * @param klabel contains superpixel labels
   */
  bool generateNeighborhoodMap(sidType* klabels,
                               const int width,
                               const int height);

  /**
   * Generate superpixel labels using the LKM library
   * @param a_img_name name of the input image
   */
  void generateSuperpixels(const char* a_img_name,
                           int superpixelStepSize = SUPERPIXEL_DEFAULT_STEP_SIZE,
                           float M = SUPERPIXEL_DEFAULT_M);

  void generateSupernodeLabels(const char* fn_annotation,
                               bool includeBoundaryLabels,
                               bool useColorImages);

  /**
   * Generate supernode labels (background, foreground, boundary) from
   * a binary image
   * @param includeBoundaryLabels=true means that the boundary class is used
   */
  void generateSupernodeLabelFromMaskImage(const char* fn_annotation,
                                           bool includeBoundaryLabels);

  void generateSupernodeLabelsFromMultiClassMaskImage(const char* fn_annotation,
                                                      map<ulong, labelType>& classIdxToLabel);

  void generateSupernodeLabelFromTextFile(const char* fn_annotation,
                                          int _nLabels);

  /**
   * Get center of the supernode with corresponding given id
   */
  bool getCenter(int supernodeId, node& center);

  sizeSliceType getWidth() { return img_width; }
  sizeSliceType getHeight() { return img_height; }

  ulong getSize() { return img_width*img_height; }

  uchar* getRawData() { return (uchar*)(img->imageData); }

  int getIntensity(int x, int y, int z = 0);

  float getAvgIntensity(sidType supernodeId);

  /**
   * Get average intensity for each channel of a given supernode
   */
  float getAvgIntensity(int supernodeId, int& r, int &g, int &b);

  IplImage* getColoredAnnotationImage(int nstates,
                                      labelType* labels,
                                      int nLabels,
                                      const map<labelType, ulong>* labelToClassIdx);

  inline ulong getIndex(int x, int y) {return (y*img_width) + x;}

  inline ulong getIndex(const node& n) {return (n.y*img_width) + n.x;}

  string getName() { return image_name; }

  int getNbChannels() { return img->nChannels; }

  ulong getNbEdges();
  ulong getNbNodes();
  ulong getNbSupernodes();

  inline probType getProb(int sid, int label, int scale = 0)
  {
    probType prob = 0;
    if(label < nLabels) {
      supernode* s = mSupernodes[sid];
      if(s->data) {
        prob = s->data->prob_estimates[label+(nLabels*scale)];
      }
    }
    return prob;
  }

  /**
   * Returns SID for given pixel (x,y)
   */
  sidType getSid(int x, int y) {
    return getSid(x,y,0);
  }


  sidType getSid(int x, int y, int z) {
    return pixelLabels[y*img_width+x];
  }

  const map<sidType, supernode* >& getSupernodes() {
    return mSupernodes;
  }

  supernode* getSupernode(sidType sid) {
    return mSupernodes[sid];
  }

  map<sidType, supernode* >* getMutableSupernodes() {
    return &mSupernodes;
  }

  virtual eSlicePType getType() { return SLICEP_SLICE; }

  /**
   * Load map containing neighborhood supernodes
   */
  bool loadNeighborhoodMap(const char* fn_neighbors);

  /**
   * Very very slow search method : should be re-implemented
   * Use getSID instead !
   */
  int search(int x, int y);

  void setImage(IplImage* _img) { img = _img; eraseImage = true; }

 private:

  bool neighborhoodMapLoaded;
  bool supernodeLabelsLoaded;
  int min_sid;
  int avgIntensity;

  /**
   * this class will free the memory allocated for the image only if
   * one of its method loaded the image
   */
  bool eraseImage;

  inline labelType computeSupernodeLabel(supernode* s, IplImage* imAnnotation);

  inline ulong computeSupernodeLabelFromMulticlassImage(sidType sid,
                                                        IplImage* mask,
                                                        map<ulong, labelType>& classIdxToLabel);


  /**
   * Generate superpixel labels using the LKM library
   * @param superpixelStepSize decides superpixel size (which will roughly be superpixelStepSize^2 pixels)
   */
  void generateSuperpixels(int superpixelStepSize = SUPERPIXEL_DEFAULT_STEP_SIZE,
                           float M = SUPERPIXEL_DEFAULT_M);

  void init(int width, int height);

  /**
   * Load superpixel labels from a given file
   * @param img_name name of the image
   * @param fn_label is a file containing the labels for each pixel
   */
  void initSuperpixels(const char* a_image_name, const char* fn_label,
                       int superpixelStepSize, float M,
                       bool bGenerateNeighborhoodMap);
};

#endif // SLICE_H
