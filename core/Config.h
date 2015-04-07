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

#ifndef CONFIG_H
#define CONFIG_H

#include <map>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

//------------------------------------------------------------------------------

enum eConfigType {
  CONFIG_FILE = 0,
  CONFIG_STRING
};

//------------------------------------------------------------------------------

class Config
{
 public:
  static Config* pInstance;

  map<string,string> parameters;

  static Config* Instance() 
  {    
    if (pInstance == 0)  // is it the first call?
      {
        printf("[Config] Error : you have to load a configuration file\n");
        exit(-1);
        //pInstance = new Config(); // create unique instance
      }
    return pInstance; // address of unique instance
  }

  Config(string input, eConfigType configType = CONFIG_FILE);

  static void setInstance(Config* aInstance) {
    pInstance = aInstance;
  }

  bool getParameter(string parameterName, string& parameterValue);

  bool readConfigFile(string filename);

  bool readConfigString(string input);

  bool setParameter(string parameterName, const string& parameterValue);
};

#endif // CONFIG_H

