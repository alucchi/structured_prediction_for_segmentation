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
// Written and (C) by Aurelien Lucchi and Kevin Smith                  //
// Contact <aurelien.lucchi@gmail.com> for comments & bug reports      //
// Contact <kevin.smith@epfl.ch> for comments & bug reports            //
/////////////////////////////////////////////////////////////////////////

#include "F_GradientStats.h"
#include "Histogram.h"
#include "utils.h"
#include "Config.h"

// standard libraries
#include <cv.h>
#include <dirent.h>
#include <highgui.h>
#include <math.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

// ITK Header Files
#include "itkImage.h"
#include "itkRecursiveGaussianImageFilter.h"
#include "itkImportImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"

using namespace std;

//-----------------------------------------------------------------------------

/**
 * @param direction : 0=X, 1=Y, 2=Z
 */
template <typename TInputPixelType, typename TOutputPixelType>
void gaussianFilter(unsigned char* rawData, long nx, long ny, long nz, long nc,
                    int direction,
                    TOutputPixelType*& output1d)
{
  typedef TInputPixelType  PixelType;

  PixelType * data1d = reinterpret_cast< PixelType * >(rawData);

  ulong pagesz = nx*ny*nz;		
  unsigned long int numberOfPixels = pagesz*nc;
  long offsets = 0;
  const unsigned int Dimension = 3;
		
  typedef itk::Image< TInputPixelType, Dimension > InputImageType;
  typedef itk::Image< TOutputPixelType, Dimension > OutputImageType;
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
  region.SetSize( size );
		
  typename InputImageType::PointType origin;
  origin.Fill( 0.0 );
		
  importFilter->SetOrigin( origin );
		
  typename ImportFilterType::SpacingType spacing;
  spacing.Fill( 1.0 );
  importFilter->SetSpacing( spacing );
		
  const bool importImageFilterWillOwnTheBuffer = false;
		
  typedef itk::CastImageFilter< InputImageType, OutputImageType> CastImageFilterType;
  typename CastImageFilterType::Pointer castImageFilter = CastImageFilterType::New();
		
  typedef itk::RecursiveGaussianImageFilter<OutputImageType, OutputImageType> GaussianImageFilterType;
  typename GaussianImageFilterType::Pointer gaussianImageFilter = GaussianImageFilterType::New();


  try
    {
      output1d = new TOutputPixelType [numberOfPixels];
    }
  catch(...)
    {
      std::cerr << "[GaussianFilter] Error while allocating memory." << std::endl;
      return;
    }
				
  const bool filterWillDeleteTheInputBuffer = false;
				
  for(long ch=0; ch<nc; ch++)
    {
      offsets = ch*pagesz;
					
      TOutputPixelType *p = output1d+offsets;
					
      importFilter->SetImportPointer( data1d+offsets, pagesz, importImageFilterWillOwnTheBuffer );
					
      castImageFilter->SetInput( importFilter->GetOutput() );
      try
        {
          castImageFilter->Update();
        }
      catch( itk::ExceptionObject & excp)
        {
          std::cerr << "Error while running this filter." << std::endl;
          std::cerr << excp << std::endl;
          return;
        }

      // set parameters
      gaussianImageFilter->SetInput( castImageFilter->GetOutput() );					
      gaussianImageFilter->SetDirection(direction);
      gaussianImageFilter->SetFirstOrder();
      gaussianImageFilter->SetSigma(3);
					
      gaussianImageFilter->GetOutput()->GetPixelContainer()->SetImportPointer( p, pagesz, filterWillDeleteTheInputBuffer);
					
      try
        {
          gaussianImageFilter->Update();
        }
      catch( itk::ExceptionObject & excp)
        {
          std::cerr << "Error while running this filter." << std::endl;
          std::cerr << excp << std::endl;
          return;
        }
					
    }
}

template <typename TInputPixelType, typename TOutputPixelType>
void gradientMagnitude(unsigned char* rawData, long nx, long ny, long nz, long nc,
                       float gaussianVariance, TOutputPixelType*& output1d)
{
  typedef TInputPixelType  PixelType;
	
  PixelType * data1d = reinterpret_cast< PixelType * >(rawData);
  ulong pagesz = nx*ny*nz;
  unsigned long int numberOfPixels = pagesz*nc;
  long offsets = 0;
  const unsigned int Dimension = 3;
		
  typedef itk::Image< TInputPixelType, Dimension > InputImageType;
  typedef itk::Image< TOutputPixelType, Dimension > OutputImageType;
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
  region.SetSize( size );
		
  typename InputImageType::PointType origin;
  origin.Fill( 0.0 );
		
  importFilter->SetOrigin( origin );
		
  typename ImportFilterType::SpacingType spacing;
  spacing.Fill( 1.0 );	
  importFilter->SetSpacing( spacing );
		
  const bool importImageFilterWillOwnTheBuffer = false;
		
  typedef itk::GradientMagnitudeRecursiveGaussianImageFilter<InputImageType, OutputImageType> GradientType;
  typename GradientType::Pointer gradientFilter = GradientType::New();

  gradientFilter->SetSigma(gaussianVariance);

  try
    {
      output1d = new TOutputPixelType [numberOfPixels];
    }
  catch(...)
    {
      std::cerr << "[F_GradientStats] Error while allocating memory." << std::endl;
      return;
    }
				
  const bool filterWillDeleteTheInputBuffer = false;
				
  for(long ch=0; ch<nc; ch++)
    {
      offsets = ch*pagesz;
      TOutputPixelType *p = output1d+offsets;
					
      importFilter->SetImportPointer( data1d+offsets, pagesz, importImageFilterWillOwnTheBuffer );
      gradientFilter->SetInput( importFilter->GetOutput() );

      gradientFilter->GetOutput()->GetPixelContainer()->SetImportPointer( p, pagesz, filterWillDeleteTheInputBuffer);
					
      try
        {
          gradientFilter->Update();
        }
      catch( itk::ExceptionObject & excp)
        {
          std::cerr << "[F_GradientStats] Error ehile running this filter." << std::endl;
          std::cerr << excp << std::endl;
          return;
        }

    }
}

//--------------------------------------------------------------------- METHODS

F_GradientStats::F_GradientStats(Slice3d& slice3d,
                                 eGradientType gradientType)
{
  /*
  float gaussianVariance = 10.0f;
  string paramGaussianVariance;
  if(Config::Instance()->getParameter("gradientStats_gaussianVariance", paramGaussianVariance)) {
    gaussianVariance = atof(paramGaussianVariance.c_str());
    PRINT_MESSAGE("[F_GradientStats] gaussianVariance=%f\n", gaussianVariance);
  }
  */


  minGradientValue = new float[numScales];
  maxGradientValue = new float[numScales];
  fGradientCube = new float*[numScales];

  nBinsPerScale = 10;
  computeGradientCube(slice3d,
                      gradientType);
}

F_GradientStats::~F_GradientStats()
{
  for(int s = 0; s < numScales; ++s) {
    delete[] fGradientCube[s];
  }
  delete[] minGradientValue;
  delete[] maxGradientValue;
  delete[] fGradientCube;
}

int F_GradientStats::getSizeFeatureVectorForOneSupernode()
{
  return nBinsPerScale*numScales;
}

bool F_GradientStats::getFeatureVectorForOneSupernode(osvm_node *x, Slice* slice, int supernodeId)
{
  // not implemented yet
  return false;
}

void F_GradientStats::computeGradientCube(Slice3d& slice3d, eGradientType gradientType)
{
  ulong cubeSize = slice3d.width*slice3d.height*(ulong)slice3d.depth;
  PRINT_MESSAGE("[F_GradientStats] Computing gradient (type = %d) for cube containing %ld voxels\n",
                (int)gradientType, cubeSize);

  for(int s = 0; s < numScales; ++s) {
    if(gradientType == GRADIENT_MAG) {
      gradientMagnitude<unsigned char, float>(slice3d.raw_data,
                                              (long)slice3d.width,(long)slice3d.height,(long)slice3d.depth,(long)slice3d.nChannels,
                                              scales[s], fGradientCube[s]);
    } else {
      gaussianFilter<unsigned char, float>(slice3d.raw_data,
                                           (long)slice3d.width,(long)slice3d.height,(long)slice3d.depth,(long)slice3d.nChannels,
                                           (int)gradientType,
                                           fGradientCube[s]);
    }

    minGradientValue[s] = fGradientCube[s][0];
    maxGradientValue[s] = fGradientCube[s][0];

#if 0
    for(sizeSliceType i = 0; i < cubeSize; i++) {
      if(minGradientValue[s] > fGradientCube[s][i])
        minGradientValue[s] = fGradientCube[s][i];
      if(maxGradientValue[s] < fGradientCube[s][i])
        maxGradientValue[s] = fGradientCube[s][i];
    }
#else
    ulong idx_min = 0.01*cubeSize;
    ulong idx_max = 0.99*cubeSize; 

    vector<float> vfGradientCube(cubeSize);
    for(sizeSliceType i = 0; i < cubeSize; i++) {
      vfGradientCube[i] = fGradientCube[s][i];
    }

    nth_element(vfGradientCube.begin(), vfGradientCube.begin()+idx_min, vfGradientCube.end());
    minGradientValue[s] = vfGradientCube[idx_min];                                                                                           
                                                                                                                                 
    nth_element(vfGradientCube.begin(), vfGradientCube.begin()+idx_max, vfGradientCube.end());
    maxGradientValue[s] = vfGradientCube[idx_max];                                                                                            
#endif

    PRINT_MESSAGE("[F_GradientStats] minGradientValue scale %f : %f\n", scales[s], minGradientValue[s]);
    PRINT_MESSAGE("[F_GradientStats] maxGradientValue scale %f : %f\n", scales[s], maxGradientValue[s]);
  }
}

bool F_GradientStats::getFeatureVectorForOneSupernode(osvm_node *x,
                                                      Slice3d* slice3d,
                                                      const int supernodeId)
{
  Histogram hist(nBinsPerScale*numScales);

  supernode* s = slice3d->getSupernode(supernodeId);
  float value;
  const ulong size_slice = slice3d->width * slice3d->height;
  node n;
  int idx;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd()) {
    ni.get(n);
    ni.next();

    for(int sc = 0; sc < numScales; ++sc) {
      double valToIdx = hist.nBins/(double)(maxGradientValue[sc]-minGradientValue[sc]);
      value = fGradientCube[sc][n.z*size_slice + n.y*slice3d->width + n.x];
      value -= minGradientValue[sc];
      idx = (int)(value*valToIdx+0.5f); // rounding
      if(idx < 0) {
        idx = 0;
      }
      if(idx >= nBinsPerScale) {
        idx = nBinsPerScale - 1;
      }
      hist.histData[(sc*nBinsPerScale) + idx]++;
    }
  }

  for(int i = 0; i < hist.nBins; i++) {
    x[i].value = hist.histData[i];
  }

  return true;
}
