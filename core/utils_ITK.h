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

#ifndef ITK_UTILS_H
#define ITK_UTILS_H

#ifdef USE_ITK

#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

// ITK Header Files
#include "itkBinaryDilateImageFilter.h"
#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h" 
#include "itkConnectedComponentFunctorImageFilter.h"
#include "itkGrayscaleErodeImageFilter.h"
#include "itkImportImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkLabelObject.h"
#include "itkLabelMap.h"
#include "itkLabelMapToBinaryImageFilter.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkRGBAPixel.h"
#include "itkShapeLabelObject.h"

// SliceMe
#include "Slice3d.h"
#include "globalsE.h"
#include "utils.h"

template <typename TInputPixelType, typename TOutputPixelType>
void connectedComponents(TInputPixelType* inputData,
                         long nx, long ny, long nz,
                         TOutputPixelType*& output1d,
                         const unsigned int minimumSize);

template <typename TInputPixelType, typename TOutputPixelType>
void connectedComponents_2d(TInputPixelType* inputData,
                            long nx, long ny, long nz,
                            TOutputPixelType*& output1d,
                            const unsigned int minimumSize);

template <typename TInputPixelType, typename LabelImageType,
          typename TOutputPixelType, unsigned int dimension>
void selectConnectedComponents(LabelImageType* labelImage,
                               TOutputPixelType* output1d,
                               const unsigned int minimumSize)
{
  typedef itk::ShapeLabelObject< long, dimension >       LabelObjectType;
  //typedef itk::ShapeLabelObject< TInputPixelType, dimension >       LabelObjectType;
  typedef itk::LabelMap< LabelObjectType >             LabelMapType;
  typedef itk::LabelImageToShapeLabelMapFilter< LabelImageType, LabelMapType > LabelFilterType;
  typedef itk::Image< TOutputPixelType, dimension > OutputImageType;
  typedef itk::LabelMapToBinaryImageFilter< LabelMapType, OutputImageType > BinarizerType;

  typename LabelFilterType::Pointer labelFilter = LabelFilterType::New();

  labelFilter->SetBackgroundValue(0);
  labelFilter->SetInput(labelImage);
  labelFilter->Update();

  LabelMapType * outputLabelMap = labelFilter->GetOutput();
  unsigned long numberOfLabelMapObjects = outputLabelMap->GetNumberOfLabelObjects();

  std::vector<LabelObjectType*> labelsToDelete;
  for( unsigned long labelCounter = 0; labelCounter < numberOfLabelMapObjects; labelCounter++ ) {
    LabelObjectType* labelObject = outputLabelMap->GetNthLabelObject( labelCounter );

    //printf("[utils_ITK] Component size=%d\n", (int)labelObject->GetSize());

    if(labelObject->GetSize() < minimumSize) {
      labelsToDelete.push_back(labelObject);
    }
    //outputLabelMap->RemoveLabelObject(labelObject);
    //outputImage->AddLabelObject(labelObject);
  }

  printf("[utils_ITK] %d/%d components to delete\n", (int)labelsToDelete.size(),numberOfLabelMapObjects);

  typename std::vector< LabelObjectType* >::iterator it;
  for(it = labelsToDelete.begin();
      it != labelsToDelete.end(); it++) {
    outputLabelMap->RemoveLabelObject(*it);
  }

  typename OutputImageType::SizeType size = labelImage->GetLargestPossibleRegion().GetSize();
  long nx = size[0];
  long ny = size[1];
  long nz = 1;
  if(dimension == 3)
    nz = size[2];
  long numberOfPixels = nx*ny*nz;

  typename BinarizerType::Pointer binarizer = BinarizerType::New();
  binarizer->SetInput(outputLabelMap);
  binarizer->SetForegroundValue(255);
  binarizer->SetBackgroundValue(0);
  binarizer->GetOutput()->GetPixelContainer()->SetImportPointer(output1d,numberOfPixels,false);
  //binarizer->SetBackgroundImage( this->GetInput() );
  binarizer->Update();

  /*
  {
    printf("[Rays3d] Writing output image\n");
    typedef itk::ImageFileWriter< OutputImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName("cannyOutput.tif");
    writer->SetInput(binarizer->GetOutput());
    writer->Update();
  }
  */
}


template <typename TInputPixelType, typename TOutputPixelType>
void connectedComponents(TInputPixelType* inputData,
                         long nx, long ny, long nz,
                         TOutputPixelType*& output1d,
                         const unsigned int minimumSize)
{
  const unsigned int dimension = 3;
  typedef itk::Image< TInputPixelType, dimension > InputImageType;
  typedef unsigned int TLabelPixelType;
  typedef itk::Image< TLabelPixelType, dimension > LabelImageType;

  typedef itk::ConnectedComponentImageFilter< InputImageType, LabelImageType > CCFilterType;
  typename CCFilterType::Pointer ccFilter;

  ccFilter = CCFilterType::New();

  // set parameters
  ccFilter->SetFullyConnected(true);
  ccFilter->SetBackgroundValue(0);

  // import data to an itk image
  typedef itk::ImportImageFilter< TInputPixelType, dimension > ImportFilterType;
		
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

  // allocate memory for output
  output1d = new TOutputPixelType [nx*ny*nz];
  //ccFilter->GetOutput()->GetPixelContainer()->SetImportPointer(output1d,0,false);

  // run filter
  ccFilter->SetInput(importFilter->GetOutput());
  ccFilter->Update();

  /*
  {
    printf("[Rays3d] Writing label image\n");
    typedef itk::ImageFileWriter< LabelImageType > WriterType;
    typename WriterType::Pointer writer = WriterType::New();
    writer->SetFileName("outputCC.mhd");
    writer->SetInput(ccFilter->GetOutput());
    writer->Update();
  }
  */

  selectConnectedComponents<TInputPixelType, LabelImageType, TOutputPixelType, 3>(ccFilter->GetOutput(),
                                                                               output1d,minimumSize);
}

template <typename TInputPixelType, typename TOutputPixelType>
void connectedComponents_2d(TInputPixelType* inputData,
                            long nx, long ny, long nz,
                            TOutputPixelType*& output1d,
                            const unsigned int minimumSize)
{
  const unsigned int dimension = 2;
  typedef itk::Image< TInputPixelType, dimension > InputImageType;
  typedef unsigned int TLabelPixelType;
  typedef itk::Image< TLabelPixelType, dimension > LabelImageType;
  typedef itk::ConnectedComponentImageFilter< InputImageType, LabelImageType > CCFilterType;
  // import data to an itk image
  typedef itk::ImportImageFilter< TInputPixelType, dimension > ImportFilterType;
		
  typename ImportFilterType::SizeType size;
  size[0] = nx;
  size[1] = ny;
  size[2] = 1;
		
  typename ImportFilterType::IndexType start;
  start.Fill( 0 );
		
  typename ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );
		
  typename InputImageType::PointType origin;
  origin.Fill( 0.0 );	       
				
  typename ImportFilterType::SpacingType spacing;
  spacing.Fill( 1.0 );

  // allocate memory for output
  output1d = new TOutputPixelType [nx*ny*nz];
  int sliceSize = nx*ny;
  for(int z=0; z < nz; z++)
    {
      typename CCFilterType::Pointer ccFilter = CCFilterType::New();
      // set parameters
      ccFilter->SetFullyConnected(true);
      ccFilter->SetBackgroundValue(0);

      typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
      importFilter->SetRegion( region );
      importFilter->SetImportPointer(inputData+(z*sliceSize), 0, false);
      importFilter->SetSpacing(spacing);
      importFilter->SetOrigin( origin );
      //ccFilter->GetOutput()->GetPixelContainer()->SetImportPointer(output1d+(z*sliceSize),0,false);

      // run filter
      ccFilter->SetInput(importFilter->GetOutput());
      ccFilter->Update();

      selectConnectedComponents<TInputPixelType, LabelImageType, TOutputPixelType, 2>(ccFilter->GetOutput(),
                                                                                      output1d+(z*sliceSize),
                                                                                      minimumSize);
    }
}

#if ITK_VERSION_MAJOR < 4
template <typename TInputPixelType>
void getPixelList_2dComponents(TInputPixelType* inputData,
                               long nx, long ny, long nz,
                               set<long>& pixelList)
{
  const unsigned int dimension = 2;
  typedef itk::Image< TInputPixelType, dimension > InputImageType;
  typedef unsigned int TLabelPixelType;
  typedef itk::Image< TLabelPixelType, dimension > LabelImageType;
  typedef itk::ConnectedComponentImageFilter< InputImageType, LabelImageType > CCFilterType;
  // import data to an itk image
  typedef itk::ImportImageFilter< TInputPixelType, dimension > ImportFilterType;
		
  typename ImportFilterType::SizeType size;
  size[0] = nx;
  size[1] = ny;
  size[2] = 1;
		
  typename ImportFilterType::IndexType start;
  start.Fill( 0 );
		
  typename ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );
		
  typename InputImageType::PointType origin;
  origin.Fill( 0.0 );	       
				
  typename ImportFilterType::SpacingType spacing;
  spacing.Fill( 1.0 );

  // allocate memory for output
  int sliceSize = nx*ny;
  for(int z=0; z < nz; z++) {
    typename CCFilterType::Pointer ccFilter = CCFilterType::New();
    // set parameters
    ccFilter->SetFullyConnected(true);
    ccFilter->SetBackgroundValue(0);

    typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
    importFilter->SetRegion( region );
    importFilter->SetImportPointer(inputData+(z*sliceSize), 0, false);
    importFilter->SetSpacing(spacing);
    importFilter->SetOrigin( origin );

    // run filter
    ccFilter->SetInput(importFilter->GetOutput());
    ccFilter->Update();

    LabelImageType* labelImage = ccFilter->GetOutput();

    typedef itk::ShapeLabelObject< long, dimension >       LabelObjectType;
    //typedef itk::ShapeLabelObject< TInputPixelType, dimension >       LabelObjectType;
    typedef itk::LabelMap< LabelObjectType >             LabelMapType;
    typedef itk::LabelImageToShapeLabelMapFilter< LabelImageType, LabelMapType > LabelFilterType;

    typename LabelFilterType::Pointer labelFilter = LabelFilterType::New();

    labelFilter->SetBackgroundValue(0);
    labelFilter->SetInput(labelImage);
    labelFilter->Update();

    LabelMapType * outputLabelMap = labelFilter->GetOutput();
    unsigned long numberOfLabelMapObjects = outputLabelMap->GetNumberOfLabelObjects();

    for( unsigned long labelCounter = 0; labelCounter < numberOfLabelMapObjects; labelCounter++ ) {
      LabelObjectType* labelObject = outputLabelMap->GetNthLabelObject( labelCounter );

      //pixelList.insert(labelObject->GetSize());
      //printf("[utils_ITK] Component size=%d\n", (int)labelObject->GetSize());

      typename LabelObjectType::LineContainerType lineContainer = labelObject->GetLineContainer();
      typename LabelObjectType::IndexType currentIdx = lineContainer.begin()->GetIndex();
      typename LabelObjectType::LengthType  currentLength = lineContainer.begin()->GetLength();

      typename LabelObjectType::LineContainerType::const_iterator it = lineContainer.begin();

      while ( it != lineContainer.end() ) {
        const typename LabelObjectType::LineType & line = *it;
        typename LabelObjectType::IndexType        idx = line.GetIndex();
        typename LabelObjectType::LengthType    length = line.GetLength();
       
        for(int i = 0; i < length; ++i) {
          pixelList.insert(idx[0]);
        }
   
        it++;
      }
    }

    /*
    // then check the lines consistancy
    // we'll proceed line index by line index
 
     */

  }
}
#endif

template <typename TInputPixelType, typename InputImageType>
  typename InputImageType::Pointer ImportFilterFromRawData(TInputPixelType* inputData,
                         long nx, long ny, long nz)
{
  const unsigned int dimension = 3;

  // import data to an itk image
  typedef itk::ImportImageFilter< TInputPixelType, dimension > ImportFilterType;		
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
  ulong totalSize = nx*ny*nz;
  importFilter->SetImportPointer(inputData, totalSize, false);

  importFilter->Update();

  return importFilter->GetOutput();
}

template <typename TInputPixelType>
void itkImportImage(const char* imgFileName,
                    TInputPixelType*& outputData,
                    int& width,
                    int& height)
{
  const int Dimension = 2;
  typedef unsigned char PixelType;
  typedef itk::Image< PixelType, Dimension >   ImageType;
  typedef itk::ImageFileReader< ImageType >  ImageReaderType;
  typedef ImageType::RegionType RegionType; 

  ImageReaderType::Pointer imageReader = ImageReaderType::New();
  imageReader->SetFileName(imgFileName);
 
  //Load it
  ImageType *img = imageReader->GetOutput();
  img->SetBufferedRegion(img->GetLargestPossibleRegion());
  try {
    imageReader->Update();
  } catch( itk::ExceptionObject & excep ) {
    cout << "[Utils] Exception Caught !" << std::endl;
    cout << excep << std::endl;
    exit(-1);
  }

  ImageType::SizeType size = img->GetLargestPossibleRegion().GetSize();
  ulong totalSize = size[0]*size[1];
  outputData = new uchar[totalSize];
  memcpy(outputData,img->GetBufferPointer(),totalSize);

  width = size[0];
  height = size[1];
}

template <typename PixelType>
void itkExportImage(const char* imgFileName,
                    PixelType* inputData,
                    int width,
                    int height)
{
  const unsigned int Dimension = 2;
  typedef itk::Image< PixelType, Dimension > ImageType;
  typedef itk::ImportImageFilter< PixelType, Dimension >   ImportFilterType;

  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
 
  typename ImportFilterType::SizeType  size;
  size[0]  = width;  // size along X
  size[1]  = height;  // size along Y
 
  typename ImportFilterType::IndexType start;
  start.Fill( 0 );
 
  typename ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );
  importFilter->SetRegion( region );
 
  double origin[ Dimension ];
  origin[0] = 0.0;    // X coordinate
  origin[1] = 0.0;    // Y coordinate
  importFilter->SetOrigin( origin );
 
  double spacing[ Dimension ];
  spacing[0] = 1.0;    // along X direction
  spacing[1] = 1.0;    // along Y direction 
  importFilter->SetSpacing( spacing );
 
  const unsigned int numberOfPixels = width*height;

  const bool importImageFilterWillOwnTheBuffer = false;
  importFilter->SetImportPointer(inputData, numberOfPixels,
                                 importImageFilterWillOwnTheBuffer );
 
  typedef itk::ImageFileWriter< ImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
 
  writer->SetFileName(imgFileName);
  writer->SetInput(  importFilter->GetOutput()  );
  writer->Update();
}

template <typename TInputPixelType, typename TOutputPixelType>
void erode(TInputPixelType* inputData,
           long nx, long ny, long nz,
           TOutputPixelType*& output1d)
{
  const unsigned int dimension = 3;
  typedef itk::Image< TInputPixelType, dimension > ImageType;
  typedef itk::ImportImageFilter< TInputPixelType, dimension > ImportFilterType;

  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
			
  typename ImportFilterType::SizeType size;
  size[0] = nx;
  size[1] = ny;
  size[2] = nz;
			
  typename ImportFilterType::IndexType start;
  start.Fill(0);
  
  typename ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );
  importFilter->SetRegion( region );                      
  region.SetSize( size );
			
  typename ImageType::PointType origin;
  origin.Fill( 0.0 );
  importFilter->SetOrigin( origin );

  typename ImportFilterType::SpacingType spacing;
  spacing.Fill( 1.0 );
  importFilter->SetSpacing( spacing );
			
  const bool importImageFilterWillOwnTheBuffer = false;
  ulong pagesz = nx*ny*nz;
  importFilter->SetImportPointer(inputData, pagesz, importImageFilterWillOwnTheBuffer );
			
  typedef itk::BinaryBallStructuringElement<TInputPixelType, dimension > StructuringElementType;
			
  typedef unsigned char CharPixelType;	
  typedef itk::Image<CharPixelType, dimension> CharImageType;	
  typedef itk::BinaryErodeImageFilter< ImageType, CharImageType, StructuringElementType > MyBinaryFilterType;
  typename MyBinaryFilterType::Pointer filter_b = MyBinaryFilterType::New();
  //typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType, StructuringElementType > MyGrayscaleFilterType;
  //typename MyGrayscaleFilterType::Pointer filter_g = MyGrayscaleFilterType::New();
  
  filter_b->SetInput( importFilter->GetOutput() );
  //filter_g->SetInput( importFilter->GetOutput() );

  //set the SE	
  StructuringElementType  structuringElement;	
  //structuringElement.SetRadius(1);  // 3x3 structuring element
  structuringElement.SetRadius(3);
  structuringElement.CreateStructuringElement();
			
  filter_b->SetKernel( structuringElement );
  //filter_g->SetKernel(  structuringElement );
			
  TInputPixelType foreground = 255;
  filter_b->SetErodeValue( foreground );

  output1d = new TOutputPixelType[pagesz];
  const bool filterWillDeleteTheInputBuffer = false;
  filter_b->GetOutput()->GetPixelContainer()->SetImportPointer(output1d, pagesz, filterWillDeleteTheInputBuffer);

  //now compute	
  try
    {
      filter_b->Update(); //the actual computation happens here!
      //filter_g->Update(); 
    }
  catch ( itk::ExceptionObject & excp )
    {
      printf("Exception!\n");
      return;
    }
}

template <typename TInputPixelType, typename TOutputPixelType>
void dilate(TInputPixelType* inputData,
            long nx, long ny, long nz,
            TOutputPixelType*& output1d)
{
  const unsigned int dimension = 3;
  typedef itk::Image< TInputPixelType, dimension > ImageType;
  typedef itk::ImportImageFilter< TInputPixelType, dimension > ImportFilterType;

  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
			
  typename ImportFilterType::SizeType size;
  size[0] = nx;
  size[1] = ny;
  size[2] = nz;
			
  typename ImportFilterType::IndexType start;
  start.Fill(0);
  
  typename ImportFilterType::RegionType region;
  region.SetIndex( start );
  region.SetSize(  size  );
  importFilter->SetRegion( region );                      
  region.SetSize( size );
			
  typename ImageType::PointType origin;
  origin.Fill( 0.0 );
  importFilter->SetOrigin( origin );

  typename ImportFilterType::SpacingType spacing;
  spacing.Fill( 1.0 );
  importFilter->SetSpacing( spacing );
			
  const bool importImageFilterWillOwnTheBuffer = false;
  ulong pagesz = nx*ny*nz;
  importFilter->SetImportPointer(inputData, pagesz, importImageFilterWillOwnTheBuffer );
			
  typedef itk::BinaryBallStructuringElement<TInputPixelType, dimension > StructuringElementType;
			
  typedef unsigned char CharPixelType;	
  typedef itk::Image<CharPixelType, dimension> CharImageType;	
  typedef itk::BinaryDilateImageFilter< ImageType, CharImageType, StructuringElementType > MyBinaryFilterType;
  typename MyBinaryFilterType::Pointer filter_b = MyBinaryFilterType::New();
  //typedef itk::GrayscaleErodeImageFilter< ImageType, ImageType, StructuringElementType > MyGrayscaleFilterType;
  //typename MyGrayscaleFilterType::Pointer filter_g = MyGrayscaleFilterType::New();
  
  filter_b->SetInput( importFilter->GetOutput() );
  //filter_g->SetInput( importFilter->GetOutput() );

  //set the SE	
  StructuringElementType  structuringElement;	
  //structuringElement.SetRadius(1);  // 3x3 structuring element
  structuringElement.SetRadius(3);
  structuringElement.CreateStructuringElement();
			
  filter_b->SetKernel( structuringElement );
  //filter_g->SetKernel(  structuringElement );
			
  TInputPixelType foreground = 255;
  filter_b->SetDilateValue( foreground );

  output1d = new TOutputPixelType[pagesz];
  const bool filterWillDeleteTheInputBuffer = false;
  filter_b->GetOutput()->GetPixelContainer()->SetImportPointer(output1d, pagesz, filterWillDeleteTheInputBuffer);

  //now compute	
  try
    {
      filter_b->Update(); //the actual computation happens here!
      //filter_g->Update(); 
    }
  catch ( itk::ExceptionObject & excp )
    {
      printf("Exception!\n");
      return;
    }
}

template <typename TInputPixelType>
void exportColorAlphaTIFCube(TInputPixelType* rawData,
                             const char* filename,
                             int cubeDepth,
                             int cubeHeight,
                             int cubeWidth)
{
  // import data to an itk image
  const int dimension = 3;
  //typedef uchar TInputPixelType;
  typedef itk::RGBAPixel<TInputPixelType> RGBAPixelType;
  typedef itk::Image< RGBAPixelType, dimension > InputImageType;
  typedef itk::Image< RGBAPixelType, dimension > OutputImageType;
  typedef itk::ImportImageFilter< RGBAPixelType, dimension > ImportFilterType;
  typename ImportFilterType::Pointer importFilter = ImportFilterType::New();
		
  typename ImportFilterType::SizeType size;
  size[0] = cubeWidth;
  size[1] = cubeHeight;
  size[2] = cubeDepth;
		
  typename ImportFilterType::IndexType start;
  start.Fill(0);
		
  typename ImportFilterType::RegionType region;
  region.SetIndex(start);
  region.SetSize(  size  );
		
  importFilter->SetRegion( region );	
		
  typename InputImageType::PointType origin;
  origin.Fill(0.0);
		
  importFilter->SetOrigin( origin );
				
  typename ImportFilterType::SpacingType spacing;
  spacing.Fill(1.0);
		
  importFilter->SetSpacing( spacing );
  RGBAPixelType* localBuffer = reinterpret_cast<RGBAPixelType* >(rawData);
  importFilter->SetImportPointer(localBuffer, 0, false);

  stringstream sout;
  sout << filename << ".nrrd";
  //printf("[Utils] Writing output cube %s\n", sout.str().c_str());
  typedef itk::ImageFileWriter< OutputImageType > WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(sout.str().c_str());
  writer->SetInput(importFilter->GetOutput());
  writer->Update();
}

#endif // USE_ITK

#endif //ITK_UTILS_H
