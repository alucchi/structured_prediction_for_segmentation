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

#ifndef STATMODEL_H
#define STATMODEL_H

// standard libraries
#include <cv.h>
#include <highgui.h>
#include <stdio.h>

//-------------------------------------------------------------------------TYPES

enum eSMType
{
    SM_HISTOGRAM = 0,
    SM_HISTOGRAM_ND,
    SM_DYNAMIC_GMM,
    SM_GMM,
    SM_COMB
};

class Slice;

//-------------------------------------------------------------------------CLASS

class StatModel
{
 public:
  eSMType type;

  virtual int predict(float* sample, double& prob)
  {
    prob = 0;
    return -1;
  }

  virtual int predict(IplImage* img, int x, int y, int imgOffset_x, int imgOffset_y, double& prob, int fdx = -1)
  {
    return predict(img, x, y, prob);
    //return -1;
  }

  virtual int predict(IplImage* img, int x, int y, double& prob, int fdx = -1)
  {
    const int numch = img->nChannels;
    float *sample = new float[numch];
    uchar* pImg = &((uchar*)(img->imageData + img->widthStep*y))[x*img->nChannels];
    for(int c = 0; c < img->nChannels; c++)
      sample[c] = pImg[c];
    delete [] sample;
    return predict(sample,prob);
  }

  /**
   * Prediction using a specified image type
   * @param imageType = {IPL_DEPTH_8U,...}
   */
  /*
  virtual int predict(IplImage* img, IplImage*& outImg, int imageType)
  {
    return -1;
  }
  */

  /**
   * Prediction
   */
  virtual int predict(IplImage* img, IplImage*& outImg, int imageType = IPL_DEPTH_8U, int fdx = -1)
  {
    double prob;
    double maxProb = 0;
    double* ptrImg;
    IplImage* outImgd = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_64F, 1);

    for(int i = 0; i < img->height; i++)
      {
        for(int j = 0; j < img->width; j++)
          {
            ptrImg = &((double*)(outImgd->imageData + outImgd->widthStep*i))[j]; //outImgd->nChannels
            predict(img,j,i,prob);
            *ptrImg = prob;

            if(maxProb < prob)
              maxProb = prob;
          }
      }

    if(imageType != IPL_DEPTH_64F)
      {
        outImg = cvCreateImage(cvSize(img->width, img->height), imageType, 1);
        //printf("[StatModel] maxProb %f\n", maxProb);
        if(fabs(maxProb) > DBL_EPSILON)
          cvConvertScale(outImgd,outImg,255.0/maxProb);
        else
          cvZero(outImg);
        cvReleaseImage(&outImgd);
      }
    else
      {
        outImg = outImgd;
      }

    return 0;
  }

  virtual int predict(IplImage* img, const char* outputFilename)
  {
    IplImage* outImg;
    predict(img,outImg);
    cvSaveImage(outputFilename, outImg);
    cvReleaseImage(&outImg);

    return 0;
  }

  /*
  virtual int predict(IplImage* img, const char* outputFilename)
  {
    double prob;
    double maxProb = 0;
    float* pf;
    IplImage* outImg = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_8U, 1);
    IplImage* outImgf = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
    //cvZero(outImgf);

    for(int i = 0; i < img->height; i++)
      {
        for(int j = 0; j < img->width; j++)
          {
            pf = &((float*)(outImgf->imageData + outImgf->widthStep*i))[j*outImgf->nChannels];
            predict(img,i,j,prob);
            pf[0] = prob;

            if(maxProb < prob)
              maxProb = prob;
          }
      }

    //cvSaveImage("DynamicGMM_output_img.png", out_img);
    printf("maxProb %f\n", maxProb);
    cvConvertScale(outImgf,outImg,255.0/maxProb);
    cvSaveImage("img.png", img);
    cvSaveImage(outputFilename, outImg);
    cvNamedWindow("outImg");
    cvShowImage("outImg",outImg);
    cvWaitKey();
    cvReleaseImage(&outImg);
    cvReleaseImage(&outImgf);
  }
  */

  /**
   * Compute average probability of a given supernode
   */
  virtual int predict(Slice* slice, int sid, int imgOffset_x, int imgOffset_y, double& prob, int fdx = -1);

  virtual int predict(Slice* slice, int sid, double& prob, int fdx = -1);

  virtual int predict(Slice* slice, const char* outputFilename, int fdx = -1);

  virtual int predict(Slice* slice, IplImage*& outImg, int fdx = -1);

  virtual void addSamples(IplImage* img,
                          IplImage* mask,
                          int maskValue)
  { return; }

  virtual void train()
  { return; }
};

#endif // STATMODEL_H
