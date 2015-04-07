
#include "Slice.h"
#include "StatModel.h"
#include "utils.h"

int StatModel::predict(Slice* slice, int sid, int imgOffset_x, int imgOffset_y, double& prob, int fdx)
{
  double pixelProb;
  prob = 0;

  supernode* s = slice->mSupernodes[sid];
  node n;
  nodeIterator ni = s->getIterator();
  ni.goToBegin();
  while(!ni.isAtEnd())
    {
      ni.get(n);
      ni.next();
      predict(slice->img, n.x,n.y,imgOffset_x,imgOffset_y, pixelProb, fdx);
      prob += pixelProb;
    }
  prob /= s->size();
  return 0;
}

int StatModel::predict(Slice* slice, int sid, double& prob, int fdx)
{
  return predict(slice, sid, 0, 0, prob, fdx);
}

int StatModel::predict(Slice* slice, const char* outputFilename, int fdx)
{
  IplImage* outImg;
  predict(slice, outImg);
  saveImage(outputFilename, outImg);
  cvReleaseImage(&outImg);
  return 0;
}

int StatModel::predict(Slice* slice, IplImage*& outImg, int fdx)
{
  double prob;
  double* ptrImg;

  // output image
  outImg = cvCreateImage(cvSize(slice->img->width, slice->img->height),
                         IPL_DEPTH_64F,1);

  supernode* s;
  node n;
  for(map<sidType, supernode* >::iterator it = slice->mSupernodes.begin();
      it != slice->mSupernodes.end(); it++)
    {
      predict(slice, it->first, prob);

      // color all nodes that belong to the supernode
      s = slice->mSupernodes[it->first];
      nodeIterator ni = s->getIterator();
      ni.goToBegin();
      while(!ni.isAtEnd())
        {
          ni.get(n);
          ni.next();

          ptrImg = &((double*)(outImg->imageData + outImg->widthStep*(n.y)))[n.x];
          *ptrImg = prob;
        }
    }
  return 0;
}
