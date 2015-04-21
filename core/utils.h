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


#include <string>
#include <stdio.h>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <dirent.h>
#include <string>
#include <cv.h>
#include <highgui.h>
#ifndef _WIN32
#include <sys/resource.h>
#endif


//SliceMe
#include "Config.h"
#include "Slice.h" // labelType
#include "Slice_P.h"
#include "globalsE.h" // sizeSliceType

// ITK Header Files
#ifdef USE_ITK
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkImportImageFilter.h"
#include "itkConnectedComponentFunctorImageFilter.h"
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkShapeLabelObject.h"
#include "itkLabelMapToBinaryImageFilter.h"
#endif

using namespace std;

#ifndef UTILS_H
#define UTILS_H

//------------------------------------------------------------------------MACROS

#ifndef MAX
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#endif

#ifdef _WIN32
#define REMOVE "del"
#define ZIP "7z"
#else
#define REMOVE "rm"
#define ZIP "zip"
#endif

//---------------------------------------------------------------------FUNCTIONS

void bresenhamLine3d(int* p1, int* p2, float*& xs, float*& ys, float*& zs,
                     int& nb_pts);

int compareBWImages(const char* imageModelName,
                    const char* imageName,
                    float& true_neg,
                    float& true_pos,
                    float& false_neg,
                    float& false_pos,
                    bool normalize = true);
/*
 * Only compare first channel of the 2 given images
 */
int compareBWImages(IplImage *ptrModel,
                    IplImage *ptrImg,
                    float& true_neg,
                    float& true_pos,
                    float& false_neg,
                    float& false_pos,
                    bool normalize = true);

void compareVolumes(uchar* annotationData,
                    uchar* data,
                    int width,
                    int height,
                    int depth,
                    float& true_neg,
                    float& true_pos,
                    float& false_neg,
                    float& false_pos,
                    bool normalize=true,
                    bool useColorAnnotations=false,
                    ulong* TP = 0,
                    ulong* TN = 0);

void compareMultiLabelVolumes(Slice_P& slice_GT,
                              const labelType* labels,
                              const int class_label,
                              float& true_neg,
                              float& true_pos,
                              float& false_neg,
                              float& false_pos,
                              bool normalize,
                              bool useColorAnnotations = false,
                              ulong* TP = 0,
                              ulong* TN = 0);

void compareMultiLabelVolumes_nodeBased(Slice_P& slice_GT,
                                        const labelType* groundtruth,
                                        const labelType* labels,
                                        const int class_label,
                                        float& true_neg,
                                        float& true_pos,
                                        float& false_neg,
                                        float& false_pos,
                                        bool normalize,
                                        bool useColorAnnotations,
                                        ulong* TP,
                                        ulong* TN);

void compareMultiLabelVolumes_givenMask_nodeBased(Slice_P& slice_GT,
                                                  const labelType* mask,
                                                  const labelType* labels,
                                                  const int class_label,
                                                  float& true_neg,
                                                  float& true_pos,
                                                  float& false_neg,
                                                  float& false_pos,
                                                  bool normalize,
                                                  bool useColorAnnotations,
                                                  ulong* TP,
                                                  ulong* TN);

bool containsImageExtension(string path);

void copyFile(const char* src_filename, const char* dst_filename);

/**
 * @param outputData memory is allocated inside the function
 */
void cubeFloat2Uchar(float* inputData, uchar*& outputData,
                     int nx, int ny, int nz);

int drawLabels(Slice* slice, const char* prediction_filename,
               const char* outputFilename, bool use_prob,
               eColorMapType colormapType);

int enumerate_files_in_dir(const char* dir,
                           vector<string> &files,
                           const char* pattern);

void exportCube(uchar* rawData,
                const char* filename,
                int cubeDepth,
                int cubeHeight,
                int cubeWidth);

void exportCube(float* rawData,
                const char* filename,
                int cubeDepth,
                int cubeHeight,
                int cubeWidth);

double getMedian(vector<double>& list_values);

bool getGroundTruthName(string& groundTruthName, const string& maskDir,
                        const string& filename);

#ifdef USE_ITK
void importTIFCube(const char* imgFileName,
                   uchar*& outputData,
                   sizeSliceType& width,
                   sizeSliceType& height,
                   sizeSliceType& depth);

void importTIFCube_noAllocation(const char* imgFileName,
                                uchar*& outputData,
                                sizeSliceType& width,
                                sizeSliceType& height,
                                sizeSliceType& depth);

void exportTIFCube(uchar* rawData,
                   const char* filename,
                   int cubeDepth,
                   int cubeHeight,
                   int cubeWidth);

void importCube(const char* imgFileName,
                uchar*& outputData,
                int& width,
                int& height,
                int& depth);

void importNRRDCube_uint(const char* imgFileName,
                         uint*& outputData,
                         int& width,
                         int& height,
                         int& depth);

void exportNRRDCube(uchar* rawData,
                    const char* filename,
                    int cubeDepth,
                    int cubeHeight,
                    int cubeWidth);

void exportNRRDCube(uint* rawData,
                    const char* filename,
                    int cubeDepth,
                    int cubeHeight,
                    int cubeWidth);

void exportColorTIFCube(uchar* rawData,
                        const char* filename,
                        int cubeDepth,
                        int cubeHeight,
                        int cubeWidth);

#endif

void exportVIVACube(float* rawData,
                    const char* filename,
                    int cubeDepth,
                    int cubeHeight,
                    int cubeWidth);

void exportVIVACube(uchar* rawData,
                    const char* filename,
                    int cubeDepth,
                    int cubeHeight,
                    int cubeWidth);

void exportImageFromCube(const char* output_name, labelType* nodeLabels,
                         int width, int height, int firstImage, int nImages);

void exportImageFromColorCube(const char* output_name, labelType* nodeLabels,
                              int width, int height, int depth, int firstImageToExport, int nImagesToExport);

bool isDirectory(const char* path);
bool isDirectory(string path);

bool fileExists(const char* filename);
bool fileExists(string filename);

string findLastFile(const string& file_pattern, const string& extension, int* _idx = 0);

string getDirectoryFromPath(string path);

string getLastDirectoryFromPath(string path);

string getNameFromPath(string path);

string getNameFromPathWithoutExtension(string path);

string getNonExistingName(string name);

/**
 * List files recursively
 */
int getFilesInDirRec(const char* dir,
                     vector<string> &files,
                     const char* ext = 0);

int getFilesInDir(const char* dir,
                  vector<string> &files,
                  const int firstIdx,
                  const char* ext = 0,
                  bool includePath = false);


int getFilesInDir(const char* dir,
                  vector<string> &files,
                  const char* ext = 0,
                  bool includePath = false);

string getExtension(string path);

void getClassToLabelMap(const char* colormapFilename,
                        map<ulong, labelType>& classIdxToLabel);

void getLabelToClassMap(const char* colormapFilename,
                        map<labelType, ulong>& labelToClassIdx);

int loadImagesInDir(const char* dir,
                    vector<IplImage*>& lImages,
                    const char* ext = 0,
                    bool includePath = false);

void save16bitsImage(const char* filename, IplImage* img);

/**
 * Load a 32 bits image to a binary file
 * Caller is responsible for freeing memory
 */
IplImage* load32bitsImage(const char* filename, CvSize& size);

/**
 * Save a 32 bits image to a binary file
 */
void save32bitsImage(const char* filename, IplImage* img);

/**
 * Load a double image to a binary file
 * Caller is responsible for freeing memory
 */
IplImage* loadDoubleImage(const char* filename, CvSize& size, int nChannels);

/**
 * Save a double image to a binary file
 */
void saveDoubleImage(const char* filename, IplImage* img);

void saveFloatImage(const char* filename, IplImage* img);

IplImage* loadFloatImage(const char* filename, CvSize& size, int nChannels);

/**
 * Save an image
 */
void saveImage(const char* filename, IplImage* img, const char* ext = ".png");

/**
 * Convert double image to uchar image
 * param imgOut memory for output image is allocated by this function
 */
void double2ucharImage(IplImage* imgIn, IplImage*& imgOut);

void float2ucharImage(IplImage* imgIn, IplImage*& imgOut);

void splitString(const string& str, vector<string>& tokens);

void splitStringUsing(const string& str, vector<string>& tokens, const char separator);

int sign(int v);

uint time_seed();

void printProcessInfo(struct rusage *p);

void crossProduct(float* a,  float* b, float* c);

float l2Norm(float* a, int n);

void matMulVec_3(float* M, float* v, float* res);

void classIdxToRGB(ulong classIdx, uchar& r, uchar& g, uchar& b);

ulong RGBToclassIdx(uchar r, uchar g, uchar b);

void zipAndDeleteCube(const char* cubeName);

void getFeatureTypes(const int featureId, vector<eFeatureType>& feature_types);

ulong getFeatureTypeId(const vector<eFeatureType>& feature_types);

void set_default_parameters(Config* config);

void getColormapName(string& paramColormap);

void loadData(string imageDir, string maskDir, Config* config, Slice_P*& slice);

void loadDataAndFeatures(string imageDir, string maskDir, Config* config, Slice_P*& slice, Feature*& feature, int* featureSize, int fileIdx = 0);

void loadFromDir(const char* dir, uchar*& raw_data,
                 int& width, int& height, int* nImgs);

void print_osvm_node(osvm_node *x, const char* title = 0);

void print_osvm_node_nz(osvm_node *x, const char* title = 0);


//------------------------------------------------------------TEMPLATE FUNCTIONS

template <typename T>
void StringCat(string* str, T t) {
  stringstream sout;
  sout << str;
  sout << t;
  *str = sout.str();
}

template <typename T>
string StringPrintf(const string& str, T t) {
  stringstream sout;
  sout << str;
  sout << t;
  return sout.str();
}

template <typename T>
string VarToString(T t) {
  stringstream sout;
  sout << t;
  return sout.str();
}

#ifdef USE_ITK


template <typename TInputPixelType, typename LabelImageType>
typename LabelImageType::Pointer getLabelImage(TInputPixelType* inputData,
                                     long nx, long ny, long nz)
{
  const unsigned int Dimension = 3;
  typedef itk::Image< TInputPixelType, Dimension > InputImageType;
  //typedef unsigned int TLabelPixelType;
  //typedef itk::Image< TLabelPixelType, Dimension > LabelImageType;

  typedef itk::ConnectedComponentImageFilter< InputImageType, LabelImageType > CCFilterType;
  typename CCFilterType::Pointer ccFilter;

  ccFilter = CCFilterType::New();

  // set parameters
  ccFilter->SetFullyConnected(true);
  ccFilter->SetBackgroundValue(0);

  // import data to an itk image
  typedef itk::ImportImageFilter< TInputPixelType, Dimension > ImportFilterType;
		
  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
		
  typename ImportFilterType::SizeType size;
  size[0] = nx;
  size[1] = ny;
  size[2] = nz;
		
  typename ImportFilterType::IndexType start;
  start.Fill( 0 );
		
  typename ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );
		
  importFilter->SetRegion( region );	
  //region.SetSize( size );
		
  typename InputImageType::PointType origin;
  origin.Fill( 0.0 );
		
  importFilter->SetOrigin( origin );
				
  typename ImportFilterType::SpacingType spacing;
  spacing.Fill( 1.0 );
		
  importFilter->SetSpacing( spacing );
  importFilter->SetImportPointer(inputData, 0, false);

  // run filter
  ccFilter->SetInput(importFilter->GetOutput());
  ccFilter->Update();

  LabelImageType* labelImage = ccFilter->GetOutput();
  return labelImage;
}

#if ITK_VERSION_MAJOR < 4
template <typename TInputPixelType>
ulong compareConnectedComponents(TInputPixelType* inputData,
                                 TInputPixelType* annotationData,
                                 long nx, long ny, long nz,
                                 ulong* nObjects)
{
  const unsigned int Dimension = 3;
  typedef unsigned int TLabelPixelType;
  typedef itk::Image< TLabelPixelType, Dimension > LabelImageType;

  LabelImageType::Pointer labelInput = getLabelImage<TInputPixelType,LabelImageType>(inputData,nx,ny,nz);
  LabelImageType::Pointer labelAnnotation = getLabelImage<TInputPixelType,LabelImageType>(annotationData,nx,ny,nz);

  typedef itk::ShapeLabelObject< long, Dimension >       LabelObjectType;
  //typedef itk::ShapeLabelObject< TInputPixelType, Dimension >       LabelObjectType;
  typedef itk::LabelMap< LabelObjectType >             LabelMapType;
  typedef itk::LabelImageToShapeLabelMapFilter< LabelImageType, LabelMapType > LabelFilterType;

  typename LabelFilterType::Pointer labelFilterInput = LabelFilterType::New();
  labelFilterInput->SetBackgroundValue(0);
  labelFilterInput->SetInput(labelInput);
  labelFilterInput->Update();

  typename LabelFilterType::Pointer labelFilterAnnotation = LabelFilterType::New();
  labelFilterAnnotation->SetBackgroundValue(0);
  labelFilterAnnotation->SetInput(labelAnnotation);
  labelFilterAnnotation->Update();

  //LabelMapType* labelMapInput = labelFilterInput->GetOutput();
  LabelMapType* labelMapAnnotation = labelFilterAnnotation->GetOutput();
  ulong nAnnotatedObjects = labelMapAnnotation->GetNumberOfLabelObjects();

  if(nObjects) {
    *nObjects = nAnnotatedObjects;
  }

  ulong sliceSize = nx*ny;

  ulong nDetectedObjects = 0;
  bool objectDetected;
  ulong cubeIndex;
  for( unsigned long labelCounter = 0; labelCounter < nAnnotatedObjects; labelCounter++ )
    {
      LabelObjectType * labelObject = labelMapAnnotation->GetNthLabelObject(labelCounter);

      objectDetected = false;

      typename LabelObjectType::LineContainerType lc = labelObject->GetLineContainer();
      typename LabelObjectType::LineContainerType::iterator ilc;
      for(ilc = lc.begin(); !objectDetected && ilc != lc.end(); ilc++)
        {
          typename LabelObjectType::LineType::IndexType idx = ilc->GetIndex();
          typename LabelObjectType::LineType::LengthType len = ilc->GetLength();
          for(int i=0;i<(int)len;i++)
            {
              cubeIndex = (idx[2]*sliceSize) + (idx[1]*nx) + (idx[0]+i);
              if(annotationData[cubeIndex])
                {
                  objectDetected = true;
                  nDetectedObjects++;
                  break;
                }
            }
        }
    }

  //printf("nDetectedObjects=%ld/%ld\n", nDetectedObjects,nAnnotatedObjects);

  return nDetectedObjects;
}
#endif

template <typename TInputPixelType>
ulong countConnectedComponents(TInputPixelType* inputData,
                               long nx, long ny, long nz)
{
  const unsigned int Dimension = 3;
  typedef itk::Image< TInputPixelType, Dimension > InputImageType;
  typedef unsigned int TLabelPixelType;
  typedef itk::Image< TLabelPixelType, Dimension > LabelImageType;

  typedef itk::ConnectedComponentImageFilter< InputImageType, LabelImageType > CCFilterType;
  typename CCFilterType::Pointer ccFilter;

  ccFilter = CCFilterType::New();

  // set parameters
  ccFilter->SetFullyConnected(true);
  ccFilter->SetBackgroundValue(0);

  // import data to an itk image
  typedef itk::ImportImageFilter< TInputPixelType, Dimension > ImportFilterType;
		
  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
		
  typename ImportFilterType::SizeType size;
  size[0] = nx;
  size[1] = ny;
  size[2] = nz;
		
  typename ImportFilterType::IndexType start;
  start.Fill( 0 );
		
  typename ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );
		
  importFilter->SetRegion( region );	
  //region.SetSize( size );
		
  typename InputImageType::PointType origin;
  origin.Fill( 0.0 );
		
  importFilter->SetOrigin( origin );
				
  typename ImportFilterType::SpacingType spacing;
  spacing.Fill( 1.0 );
		
  importFilter->SetSpacing( spacing );
  importFilter->SetImportPointer(inputData, 0, false);

  // run filter
  ccFilter->SetInput(importFilter->GetOutput());
  ccFilter->Update();

  return (ulong) ccFilter->GetObjectCount ();

  /*
  LabelImageType* labelImage = ccFilter->GetOutput();

  typedef itk::ShapeLabelObject< long, Dimension >       LabelObjectType;
  //typedef itk::ShapeLabelObject< TInputPixelType, Dimension >       LabelObjectType;
  typedef itk::LabelMap< LabelObjectType >             LabelMapType;
  typedef itk::LabelImageToShapeLabelMapFilter< LabelImageType, LabelMapType > LabelFilterType;

  typename LabelFilterType::Pointer labelFilter = LabelFilterType::New();
  labelFilter->SetBackgroundValue(0);
  labelFilter->SetInput(labelImage);
  labelFilter->Update();

  LabelMapType * outputLabelMap = labelFilter->GetOutput();
  unsigned long numberOfLabelMapObjects = outputLabelMap->GetNumberOfLabelObjects();
  return numberOfLabelMapObjects;
*/
}

#endif // USE_ITK

#endif //UTILS_H
