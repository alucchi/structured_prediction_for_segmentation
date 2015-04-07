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

#include "Histogram.h"

// std libraries
#include <fstream>

using namespace std;

Histogram::Histogram(int _nBins)
{
  init(_nBins);
}

void Histogram::init(int _nBins)
{
  nBins = _nBins;
  histData = new double[nBins];	    
  for (int i = 0; i < nBins; i++)
    histData[i] = 0;
}

void Histogram::load(const char* filename)
{
  ifstream ifs(filename);
  
  const int MAX_LENGTH = 100;
  char line[MAX_LENGTH];
  ifs.getline(line, MAX_LENGTH); // first line contains number of bins
  int _nBins = atoi(line);
  init(_nBins);

  // load bin values
  int idx = 0;
  while(ifs.getline(line, MAX_LENGTH)) {
    histData[idx] = atof(line);
    idx++;
  }

  ifs.close();
}

void Histogram::save(const char* filename)
{
  ofstream ofs(filename);
  
  ofs << nBins << endl;

  // load bin values
  for (int x = 0; x < nBins; x++) {
    ofs << histData[x] << endl;
  }

  ofs.close();
}

void Histogram::normalize() 
{
  double nEntries = getNbEntries();
  if(nEntries>1e-6)
    for (int x = 0; x < nBins; x++) {
      histData[x] = histData[x] / nEntries;
    }
}

double Histogram::getValue(int bin) 
{
  return histData[bin];
}

Histogram::~Histogram()
{
  delete[] histData;
}

void Histogram::reset()
{
  for (int i = 0; i < nBins; i++)
    histData[i] = 0;
}

double Histogram::getNbEntries()
{
  double n =0;
  for (int i = 0; i < nBins; i++)
    n += histData[i];
  return n;
}

void Histogram::print() 
{
  for (int x = 0; x < nBins; x++) {
    printf("%d:%.2g ", x,histData[x]);
  }
  printf("\n");
}

void Histogram::print(int nBinsPerLine) 
{
  printf("\t");
  for(int i=0; i < nBinsPerLine; i++)
    printf("%d\t", i);
  //printf("\n");

  int i = 0;
  for (int x = 0; x < nBins; x++) {
    if(x%nBinsPerLine==0)
      printf("\n%d",i++);
    printf("\t%.2g", histData[x]);
  }
  printf("\n");

  /*
  for (int x = 0; x < nBins; x++) {
    printf("%d:%.2g ", x,histData[x]);
    if((x+1)%nBinsPerLine==0)
      printf("\n");
  }
  printf("\n");
  */
}

void Histogram::draw(IplImage* histImage)
{
  //set all histogram values to 255
  cvSet( histImage, cvScalarAll(255), 0 );
  //create a factor for scaling along the width
  int bin_w = cvRound((double)histImage->width/nBins);

  for(int i = 0; i < nBins; i++) {
    //draw the histogram data onto the histogram image
    cvRectangle( histImage, cvPoint(i*bin_w, histImage->height),
                 cvPoint((i+1)*bin_w,
                         histImage->height - histData[i]),
                 cvScalarAll(0), -1, 8, 0 );
  }
}

