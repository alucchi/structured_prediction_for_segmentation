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

#ifndef SLICE3D_H
#define SLICE3D_H

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cv.h>
#include <highgui.h>

// SliceMe
#include "globalsE.h"
#include "Slice.h"
#include "Slice_P.h"
#include "utils.h"

using namespace std;

//#define USE_REVERSE_INDEXING // reverse indexing is memory consumming

//#define UNITIALIZED_SIZE 0
#define UNITIALIZED_SIZE -1

//--------------------------------------------------------------------- CLASSES


/*
 * Slice3d is a class used to store supervoxels
 * Each supernode is composed of nodes (pixels).
 */
class Slice3d : public Slice_P
{
 public:
  sizeSliceType width;
  sizeSliceType height;
  sizeSliceType depth;
  int numlabels;
  int nChannels;
  sizeSliceType sliceSize;
  int maxDegree; // maximum degree in the graph

#ifdef USE_REVERSE_INDEXING
  sidType** klabels;
#endif

  /**
   * Volume data
   * Order : z,y,x
   */
  uchar* raw_data;
  
  // map containing a list of supernodes indexed by their supernode id
  map<sidType, supernode* >* mSupervoxels;

  /**
   * No initialization. Should be used when importing data with the importData method.
   */
  Slice3d();

  Slice3d(uchar* a_raw_data,
          int awidth, int aheight,
          int adepth,
          sizeSliceType vstep,
          int anChannels = 1,
          bool _loadNeighbors = true);

  Slice3d(const char* input_dir,
          int vstep = DEFAULT_VOXEL_STEP,
          int nImgs = -1,
          bool _loadNeighbors = true);

  Slice3d(const char* input_dir,
          int awidth,
          int aheight,
          int adepth,
          int vstep = DEFAULT_VOXEL_STEP,
          bool _loadNeighbors = true);

  Slice3d(const char* input_dir,
          const node& _start,
          const node& _end,
          int vstep = DEFAULT_VOXEL_STEP,
          bool _loadNeighbors = true);

  ~Slice3d();

  /**
   * Return pixel intensity at the given coordinate (x,y,z)
   */
  uchar at(int x, int y, int z);

  /**
   * Compute the label of a given supernode
   * @param mask_data cubes containing the supervoxel labels
   */
  inline labelType computeSupernodeLabel(sidType sid, uchar* mask_data);

  /**
   * Compute the label of a given supernode.  This function should be
   * used when colored training images are provided.
   * @param mask_data cubes containing the supervoxel labels
   */
  labelType computeSupernodeLabel_advanced(sidType sid, uchar* mask_data);

  /**
   * Create volume with colored node labels
   * Memory is allocated by this function but caller is then responsible
   * for freeing the memory.
   */
  uchar* createNodeLabelVolume();

  /**
   * Create indexing structures named mSupervoxels
   */
  void createIndexingStructures(sidType** _klabels, bool force = false);

  void createOverlayAnnotationImage(const char* filename, int imageId);

  void createReverseIndexing(sidType**& _klabels);

  /**
   * Caller is responsible for freeing memory
   */
  void createSupernodeLabels(const uchar* nodeLabels, labelType*& labelCube,
                             int nClasses);

  bool importData(const char* filename,
                  const int iwidth = -1,
                  const int iheight = -1,
                  const int idepth = -1);

  /**
   * Import supervoxels from a text file with one label per line
   * Assumes that the correct width, height and depth were passed
   * to the constructor
   * importSupervoxelsFromBinaryFile is more efficient...
   */
  void importSupervoxels(const char* filename);

  /**
   * Import supervoxels from a binary file
   * Assumes that the correct width, height and depth were passed
   * to the constructor
   */
  void importSupervoxelsFromBinaryFile(const char* filename);

  void importSupervoxelsFromBuffer(const uint* buffer, int _width, int _height, int _depth);

  void init();

  /**
   * Export data
   */
  void exportData(const char* filename);

  void exportOverlay(const char* filename);

  void exportOverlay(const char* filename, labelType* labels);

  void exportProbabilities(const char* filename, int nClasses,
                           float* pbs);

  /**
   * Export supervoxels
   */
  void exportSupervoxels(const char* filename);

  /**
   * Export supervoxels to a binary file
   */
  void exportSupervoxelsToBinaryFile(const char* filename);

  void exportSupervoxelsToNRRD(const char* filename);

  /**
   * Export labels.
   * See exportProbabilities to export probabilities associated to each label.
   */
  void exportSupernodeLabels(const char* filename, int nClasses = 2);

  void exportSupernodeLabels(const char* filename, int nClasses,
                             vector<labelType>& labels,
                             int nLabels = -1);

  void exportSupernodeLabels(const char* filename, int nClasses,
                             labelType* labels,
                             int nLabels,
			     const map<labelType, ulong>* labelToClassIdx);

  void generateSupernodeLabels(const char* fn_annotation,
                               bool includeBoundaryLabels,
                               bool useColorImages);

  /**
   * Generate supernode labels (background, foreground, boundary) from
   * a binary ground truth cube
   * @param includeBoundaryLabels=true means that the boundary class is used
   */
  void generateSupernodeLabelFromMaskDirectory(const char* mask_dir,
                                               bool includeBoundaryLabels,
                                               bool useColorImages = false);

  /**
   * Generate supernode labels (background, foreground, boundary) from
   * a binary ground truth cube
   * @param includeBoundaryLabels=true means that the boundary class is used
   */
  void generateSupernodeLabelFromMaskImages(uchar* mask_data,
                                            bool includeBoundaryLabels,
                                            bool useColorImages = false);

  string getName() { return getLastDirectoryFromPath(inputDir); }

  inline ulong getIndex(int x, int y, int z) {return z*sliceSize+y*width+x;}

  inline ulong getIndex(const node& n) {return (n.z*sliceSize) + (n.y*width) + n.x;}

  /**
   * This code does not support color cubes.
   */
  int getNbChannels() { return 1; }

  uchar* getRawData() { return raw_data; }

#ifdef USE_REVERSE_INDEXING
  sidType getSid(int x, int y, int z) { return klabels[z][(y*width) + x]; }
#else
  sidType getSid(int x, int y, int z) { printf("[Slice3d] USE_REVERSE_INDEXING not defined\n"); assert(0); return 0; }
#endif

  ulong getSize() { return sliceSize*depth; }

  uchar* getSupernodeProbabilities(int class_nr);

  labelType* getSupernodeLabels();

  const map<sidType, supernode* >& getSupernodes() {
    return *mSupervoxels;
  }

  map<sidType, supernode* >* getMutableSupernodes() {
    return mSupervoxels;
  }

  /**
   * Load supervoxels using default parameters
   */
  void loadSupervoxels(const char* imageDir);

  void loadSupervoxels(const char* imageDir, const int voxel_step, const float _cubeness);

  void setMutableSupernodes(map<sidType, supernode* >* _supernodes) {
    mSupervoxels = _supernodes;
  }

  /**
   * conversion from raw (uchar, 1 channel) to double (1 channel)
   */
  int raw2Double(double**& ptr_data);

  /**
   * Supervoxel library needs a cube made of ints so we have to convert the cube
   * Ask for enough memory for the texels and make sure we got it before proceeding
   */
  int raw2RGB(unsigned int**& ptr_data);

  void resize(sizeSliceType w, sizeSliceType h, sizeSliceType d,
              map<sidType, sidType>* sid_mapping);

  void generateSupervoxels(const double _cubeness = 20);

  int getIntensity(int x, int y, int z = 0);

  float getAvgIntensity(sidType supernodeId);

  float getAvgIntensity(int supernodeId, int& r, int &g, int &b) {
    return -1;
  }

  sizeSliceType getWidth() { return width; }
  sizeSliceType getHeight() { return height; }
  sizeSliceType getDepth() { return depth; }

  inline probType getProb(int sid, int label, int scale = 0) {
    probType prob = 0;
    if(label < nLabels) {
      supernode* s = (*mSupervoxels)[sid];
      if(s->data) {
        prob = s->data->prob_estimates[label+(nLabels*scale)];
      }
    }
    return prob;
  }

  supernode* getSupernode(sidType sid);

  bool isSupernodeLabelsLoaded() { return supernodeLabelsLoaded; }

  void loadFromDir(const char* input_dir, int nImgs = -1);

  void loadFromDir(const char* dir, const node& start, const node& end);

  void loadFromDir(const char* dir, uchar*& raw_data,
                   int& width, int& height, int* nImgs);

  ulong getNbEdges() { return nbEdges; }

  ulong getNbNodes() { return width*height*depth; }

  ulong getNbSupernodes() { return mSupervoxels->size(); }

  void setVoxelStep(int _voxel_step) { supernode_step = _voxel_step; }

  bool getIncludeOtherLabel() { return includeOtherLabel; }

  void setIncludeOtherLabel(bool _value) { includeOtherLabel = _value; }

  float getMinPercentToAssignLabel() { return minPercentToAssignLabel; }
  void setMinPercentToAssignLabel(float _minPercentToAssignLabel) { minPercentToAssignLabel = _minPercentToAssignLabel; }

  virtual eSlicePType getType() { return SLICEP_SLICE3D; }

  void rescaleRawData();

  void setDeleteRawData(bool _val) { delete_raw_data = _val; }

  // no need to free memory as functions in Supernode will not reallocate memory.
  void unloadSupernodeLabels() { supernodeLabelsLoaded = false; }

#ifdef USE_REVERSE_INDEXING

  sidType getSID(uint x,uint y,uint z);

#endif


 private:

  bool supernodeLabelsLoaded;
  bool loadNeighbors;
  bool delete_raw_data;
  int start_x;
  int start_y;
  int start_z;

};

#endif // SLICE3D_H
