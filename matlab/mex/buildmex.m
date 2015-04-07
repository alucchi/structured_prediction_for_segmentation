
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change this
OPENCV_PATH='C:\Users\monso\code\aurelien\superpixels\branches\porting\lib\opencv-2.4.5';
OPENCV_LIB_DIR= [OPENCV_PATH '\build\lib\Release'];
%C:\Users\monso\code\aurelien\superpixels\branches\porting\lib\opencv-2.4.5\modules\core\include\opencv2\core
ITK_PATH='C:\Users\monso\code\InsightToolkit-3.20.1';

OPENCV_PATH='/home/lucchi/src/ThirdParty/OpenCV-2.4.2/';
OPENCV_LIB_DIR= [OPENCV_PATH 'lib/'];
ITK_PATH='/home/lucchi/src/ThirdParty/InsightToolkit-3.20.1/';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% includes
ROOT_PATH=[pwd '/../../'];
if ispc
    SLIC_PATH=[ROOT_PATH 'lib/slic/build/Release'];
    SLIC_INCLUDE_DIR=[ROOT_PATH 'lib/slic/'];
    SLICEME_PATH=[ROOT_PATH 'build/Release/'];
    SLICEME_INCLUDE_DIR=[ROOT_PATH 'core/'];
    DAI_PATH = [ROOT_PATH 'lib/libDAI024/build/Release/'];
    ZLIB_LIB_PATH = [ROOT_PATH 'lib/zlib-1.2.3-lib/lib/'];
    ZLIB_INCLUDE_PATH = [ROOT_PATH 'lib/zlib-1.2.3-lib/include/'];
else
    SLIC_PATH=[ROOT_PATH 'lib/slic/'];
    SLIC_INCLUDE_DIR=[ROOT_PATH 'lib/slic/'];
    SLICEME_PATH=[ROOT_PATH 'core/'];
    SLICEME_INCLUDE_DIR=[ROOT_PATH 'core/'];
    DAI_PATH = [ROOT_PATH 'lib/libDAI-0.2.4/lib/'];
    ZLIB_LIB_PATH = '';
    ZLIB_INCLUDE_PATH = '';
end

INCLUDE_PATH=['-I' SLIC_INCLUDE_DIR ' -I' SLICEME_INCLUDE_DIR ' -I' OPENCV_PATH '/include/ -I' ...
              OPENCV_PATH '/include/opencv/ -I' OPENCV_PATH '/include/opencv2/'];

if ispc
    INCLUDE_PATH=[INCLUDE_PATH ' -I' ...
                  OPENCV_PATH '\modules\core\include\ -I' OPENCV_PATH '\modules\video\include\ -I' ...
                  OPENCV_PATH '\modules\flann\include\ -I' OPENCV_PATH '\modules\calib3d\include\ -I' ...
                  OPENCV_PATH '\modules\objdetect\include\ -I' OPENCV_PATH '\modules\legacy\include\ -I' ...
                  OPENCV_PATH '\modules\highgui\include\ -I' ...
                  OPENCV_PATH '\modules\imgproc\include -I' OPENCV_PATH '\modules\features2d\include\'];
end

INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/Review'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/gdcm/src'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build/Utilities/gdcm'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build/Utilities/vxl/core'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build/Utilities/vxl/vcl'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build/Utilities/vxl/v3p/netlib'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/vxl/core'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/vxl/vcl'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/vxl/v3p/netlib'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/Numerics/Statistics'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build/Utilities'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/itkExtHdrs'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/nifti/znzlib'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/nifti/niftilib'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/expat'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build/Utilities/expat'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build/Utilities/DICOMParser'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/DICOMParser'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build/Utilities/NrrdIO'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/NrrdIO'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Utilities/MetaIO'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/SpatialObject'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/Numerics/NeuralNetworks'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/Numerics/FEM'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/IO'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/Numerics'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/Common'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/BasicFilters'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/Code/Algorithms'];
INCLUDE_PATH=[INCLUDE_PATH ' -I' ITK_PATH '/build'];

% libraries
if isunix
   LIBS=['-L../../ -lsliceme -L ' OPENCV_PATH '/lib/ -lopencv_core -lopencv_imgproc -lopencv_highgui -lopencv_ml -lopencv_video -lopencv_features2d -lopencv_calib3d -lopencv_objdetect -lopencv_contrib -lopencv_legacy -lopencv_flann -L' SLIC_PATH ' -l supervoxel  -L../../lib/libDAI-0.2.4/lib/ -ldai -lz'];
   LIBS=[LIBS ' -L' ITK_PATH '/build/bin -lITKIO -lITKStatistics -lITKNrrdIO -litkgdcm -litkjpeg12 -litkjpeg16 -litkopenjpeg -litkpng -litktiff -litkjpeg8 -lITKSpatialObject -lITKMetaIO -lITKDICOMParser -lITKEXPAT -lITKniftiio -lITKznz -litkzlib -lITKCommon -litksys -litkvnl_inst -litkvnl_algo -litkvnl -litkvcl -litkv3p_lsqr -lpthread -lm -ldl -litkNetlibSlatec -litkv3p_netlib'];
else
   LIBS=['-L../../ -L' SLICEME_PATH ' -lsliceme -L' OPENCV_LIB_DIR ' -lopencv_core245 -lopencv_imgproc245 ' ...
    '-lopencv_highgui245 -lopencv_ml245 -lopencv_video245 -lopencv_features2d245 -lopencv_calib3d245 ' ...
    '-lopencv_objdetect245 -lopencv_contrib245 -lopencv_legacy245 -lopencv_flann245 -L' ...
    SLIC_PATH ' -lsupervoxel  -L' DAI_PATH ' -ldai -L' ZLIB_LIB_PATH ' -lzlib'];
   LIBS=[LIBS ' -L' ITK_PATH '/build/bin/ -L' ITK_PATH '/build/bin/Release/ -lITKIO ' ...
    '-lITKStatistics -lITKNrrdIO -litkgdcm -litkjpeg12 -litkjpeg16 -litkopenjpeg ' ... 
    '-litkpng -litktiff -litkjpeg8 -lITKSpatialObject -lITKMetaIO -lITKDICOMParser ' ...
    '-lITKEXPAT -lITKniftiio -lITKznz -litkzlib -lITKCommon -litksys -litkvnl_inst ' ...
    '-litkvnl_algo -litkvnl -litkvcl -litkv3p_lsqr -lm -ldl -litkNetlibSlatec ' ...
    '-litkv3p_netlib'];
end

% add openmp flag
if isunix
    %FLAGS=['CFLAGS="$CFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp"'];
    FLAGS=['CFLAGS="\$CFLAGS -fopenmp" LDFLAGS="\$LDFLAGS -fopenmp"'];
elseif ispc
    FLAGS=['COMPFLAGS="$COMPFLAGS /openmp  /DWIN32 /D_WINDOWS /W3 /Zm1000 /GR /EHsc" LINKFLAGS="$LINKFLAGS  /STACK:10000000 /machine:x64 "'];
else
    warning('Mac not supported yet')
end

% compile
cmd = ['mex -v ' FLAGS ' ' INCLUDE_PATH ' ' LIBS ' mex_predict.cpp']
eval(cmd);
