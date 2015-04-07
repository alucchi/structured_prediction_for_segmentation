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

#include "F_Dft.h"

using namespace cv;

F_Dft::F_Dft(const char* image_name,
             int _n_rows, int _n_cols)
{
  n_cols = _n_cols;
  n_rows = _n_rows;

  Mat img = imread(image_name, CV_LOAD_IMAGE_GRAYSCALE);
  if( img.empty() ) {
    printf("Cannot read image file: %s\n", image_name);
    return;
  }
  int M = getOptimalDFTSize( img.rows );
  int N = getOptimalDFTSize( img.cols );
  //printf("M %d N %d\n", M, N);
  Mat padded;
  copyMakeBorder(img, padded, 0, M - img.rows, 0, N - img.cols, BORDER_CONSTANT, Scalar::all(0));

  Mat planes[] = {Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F)};
  Mat complexImg;
  merge(planes, 2, complexImg);

  dft(complexImg, complexImg);

  // compute log(1 + sqrt(Re(DFT(img))**2 + Im(DFT(img))**2))
  split(complexImg, planes);
  magnitude(planes[0], planes[1], planes[0]);
  planes[0].copyTo(mag);
  mag += Scalar::all(1);
  log(mag, mag);

  // crop the spectrum, if it has an odd number of rows or columns
  //mag = mag(Rect(0, 0, mag.cols & -2, mag.rows & -2));
  mag = mag(Rect(0, 0, n_rows, n_cols));

  /*
  int cx = mag.cols/2;
  int cy = mag.rows/2;

  // rearrange the quadrants of Fourier image
  // so that the origin is at the image center
  Mat tmp;
  Mat q0(mag, Rect(0, 0, cx, cy));
  Mat q1(mag, Rect(cx, 0, cx, cy));
  Mat q2(mag, Rect(0, cy, cx, cy));
  Mat q3(mag, Rect(cx, cy, cx, cy));

  q0.copyTo(tmp);
  q3.copyTo(q0);
  tmp.copyTo(q3);

  q1.copyTo(tmp);
  q2.copyTo(q1);
  tmp.copyTo(q2);
  */

  normalize(mag, mag, 0, 1, CV_MINMAX);
}

F_Dft::~F_Dft()
{
}

bool F_Dft::getFeatureVectorForOneSupernode(osvm_node *n, Slice* slice, int supernodeId)
{
  /*
  int featId = 0;
  for(int i = 0; i < n_rows; ++i) {
    for(int j = 0; j < n_cols; ++j) {
      n[featId].value = mag.at<double>(j,i);
      ++featId;
      //printf("%d ", featId);
    }
  }
  */

  double mean = 0;
  for(int i = 0; i < n_rows; ++i) {
    for(int j = 0; j < n_cols; ++j) {
      mean += mag.at<double>(j,i);
    }
  }
  mean /= n_cols*n_rows;

  double variance = 0;
  for(int i = 0; i < n_rows; ++i) {
    for(int j = 0; j < n_cols; ++j) {
      double t = mag.at<double>(j,i) - mean;
      variance += t*t;
    }
  }
  variance /= n_cols*n_rows;

  n[0].value = mean;
  n[1].value = variance;

  return true;
}

int F_Dft::getSizeFeatureVectorForOneSupernode()
{
  //return n_rows*n_cols;
  return 2;
}
