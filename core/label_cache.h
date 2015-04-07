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

#ifndef LABEL_CACHE_H
#define LABEL_CACHE_H

#include "svm_struct_api_types.h"

//------------------------------------------------------------------------------

class LabelCache
{
 public:
  static LabelCache* pInstance;

  static LabelCache* Instance() 
  {    
    if (pInstance == 0)  // is it the first call?
      {
        pInstance = new LabelCache; // create unique instance
      }
    return pInstance; // address of unique instance
  }

  static void setInstance(LabelCache* aInstance) {
    pInstance = aInstance;
  }

  ~LabelCache();

  bool exists(int id);

  bool getLabel(int id, LABEL& l);

  void clear();

  void setLabel(int id, LABEL& l);

 private:
  map<int, LABEL> labels;

};

#endif //LABEL_CACHE_H
