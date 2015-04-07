%function L = readSupervoxelCube(labelFilenm, s)
%% readSupervoxelCube reads supervoxels data files
%   
%   L = readSupervoxelCube(labelFilenm, size) reads a supervoxel label data file 
%   and converts it into a matlab label matrix. Size of the image
%   must be specified in SIZE.
%

%   Copyright © 2009 Computer Vision Lab, 
%   École Polytechnique Fédérale de Lausanne (EPFL), Switzerland.
%   All rights reserved.
%
%   Authors:    Kevin Smith         http://cvlab.epfl.ch/~ksmith/
%               Aurelien Lucchi     http://cvlab.epfl.ch/~lucchi/
%
%   This program is free software; you can redistribute it and/or modify it 
%   under the terms of the GNU General Public License version 2 (or higher) 
%   as published by the Free Software Foundation.
%                                                                     
% 	This program is distributed WITHOUT ANY WARRANTY; without even the 
%   implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
%   PURPOSE.  See the GNU General Public License for more details.

labelFilenm = 'cube1';

if 1
   fid = fopen([labelFilenm '.nfo']);
   tline = fgets(fid);
   while ischar(tline)
      % disp(tline)
      tline = fgets(fid);
      if(strncmp(tline,'cubeDepth',9) == 1)
         cubeDepth = str2double(tline(11:length(tline)))
      else if(strncmp(tline,'cubeHeight',9) == 1)
         cubeHeight = str2double(tline(11:length(tline)))
      else if(strncmp(tline,'cubeWidth',9) == 1)
         cubeWidth = str2double(tline(11:length(tline)))
      end
      end
      end
   end
   fclose(fid);
   s = [cubeWidth cubeHeight];
else
   s = [584 388];
   cubeDepth = 8;
end

fid = fopen(labelFilenm,'r');

for i=1:cubeDepth
    L{i} = fread(fid,[s(1) s(2)],'int');
    L{i} = double(L{i});
    %imagesc(L{i});
    %pause;
end

% Show first slice
imagesc(L{1});

fclose(fid);
