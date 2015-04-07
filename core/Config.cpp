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

#include "Config.h"
#include <fstream>
#include <string.h>
#include <vector>

#include "utils.h"

Config* Config::pInstance = 0; // initialize pointer

bool Config::readConfigFile(string filename)
{
  ifstream ifs(filename.c_str(), ios::in);
  if(!ifs) {
    printf("[Config] Error while loading %s\n",filename.c_str());
    return false;
  }

  const int MAX_LENGTH = 1000;
  char line[MAX_LENGTH];

  parameters.clear();
  char *ptr;
  while(ifs.getline(line, MAX_LENGTH)) {
    ptr=strchr(line,'=');
    if(!ptr) {
      break;
    }

    *ptr = 0;
    string key(line);
    string value(ptr+1);
    
    parameters[key] = value;
  }
  
  ifs.close();
  return true;
}

bool Config::readConfigString(string input)
{
  const char line_separator = '\n';
  const char keyvalue_separator = '=';
  vector<string> tokens;
  splitStringUsing(input, tokens, line_separator);

  parameters.clear();
  for(vector<string>::iterator it = tokens.begin(); it != tokens.end(); ++it) {
    vector<string> key_value;
    splitStringUsing(*it, key_value, keyvalue_separator);
    parameters[key_value[0]] = key_value[1];
  }

  return true;
}

Config::Config(string input, eConfigType configType)
{
  if(configType == CONFIG_FILE) {
    printf("[Config] filename = %s\n", input.c_str());
    readConfigFile(input);
  } else {
    printf("[Config] Reading string\n");
    readConfigString(input);
  }
}

bool Config::getParameter(string parameterName, string& parameterValue)
{
  bool retValue = false;
  map<string, string>::iterator iParam = parameters.find(parameterName);
  if(iParam != parameters.end()) {
    parameterValue = iParam->second;
    retValue = true;
  } else {
    parameterValue = "";
  }
  return retValue;
}

bool Config::setParameter(string parameterName, const string& parameterValue)
{
  bool retValue = false;
  map<string, string>::iterator iParam = parameters.find(parameterName);
  if(iParam != parameters.end()) {
    parameters[parameterName] = parameterValue;
    retValue = true;
  }
  return retValue;
}
