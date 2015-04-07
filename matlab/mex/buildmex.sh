#!/bin/bash

################################################################################

# Change this
OPENCV_PATH='/home/lucchi/src/ThirdParty/OpenCV-2.4.2/'
ITK_PATH="/home/lucchi/src/ThirdParty/InsightToolkit-3.20.1/"

################################################################################

if [ $# -gt 0 ]
then
    if [ $1 == "clean" ]
    then
        echo 'Cleaning project'
        make clean
        rm *.o
	    rm bin/*.mexglx
    fi
fi

# Set the include path for mex
if [ -d /usr/local/extern/include/ ]
then
    MEXPATH='/usr/local/extern/include/'
elif [ -d /usr/local/matlab/extern/include/ ]
then
    MEXPATH='/usr/local/matlab/extern/include/'
elif [ -d /usr/local/MATLAB/R2012a/extern/include/ ]
then
    MEXPATH='/usr/local/MATLAB/R2012a/extern/include/'
else
    MEXPATH='/usr/local/MATLAB/R2013a/extern/include/'
    #echo 'Error : please set the MEXPATH variable manually'
    #exit -1
fi

# includes
SLIC_PATH='../../lib/slic/'
SLICEME_PATH='../../core/'
INCLUDE_PATH="-I${SLIC_PATH} -I${SLICEME_PATH} -I../../lib/libsvm-3.0/ -I${OPENCV_PATH}/include/ -I${OPENCV_PATH}/include/opencv/"

INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/Review"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/gdcm/src"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build/Utilities/gdcm"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build/Utilities/vxl/core"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build/Utilities/vxl/vcl"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build/Utilities/vxl/v3p/netlib"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/vxl/core"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/vxl/vcl"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/vxl/v3p/netlib"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/Numerics/Statistics"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build/Utilities"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/itkExtHdrs"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/nifti/znzlib"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/nifti/niftilib"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/expat"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build/Utilities/expat"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build/Utilities/DICOMParser"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/DICOMParser"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build/Utilities/NrrdIO"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/NrrdIO"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Utilities/MetaIO"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/SpatialObject"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/Numerics/NeuralNetworks"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/Numerics/FEM"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/IO"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/Numerics"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/Common"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/BasicFilters"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/Code/Algorithms"
INCLUDE_PATH="${INCLUDE_PATH} -I ${ITK_PATH}/build"

# libraries

LIBS="-L../../ -lsliceme -L ${OPENCV_PATH}/lib/ `pkg-config --libs opencv` -L${SLIC_PATH} -l supervoxel  -L../../lib/libDAI-0.2.4/lib/ -ldai -lz"
LIBS="${LIBS} -L${ITK_PATH}/build/bin -lITKIO -lITKStatistics -lITKNrrdIO -litkgdcm -litkjpeg12 -litkjpeg16 -litkopenjpeg -litkpng -litktiff -litkjpeg8 -lITKSpatialObject -lITKMetaIO -lITKDICOMParser -lITKEXPAT -lITKniftiio -lITKznz -litkzlib -lITKCommon -litksys -litkvnl_inst -litkvnl_algo -litkvnl -litkvcl -litkv3p_lsqr -lpthread -lm -ldl -litkNetlibSlatec -litkv3p_netlib"

#GCC=/usr/bin/c++
GCC=g++
MEX_ARG="-cxx -O"
MEX_EXE=`which mex`
if [ "$MEX_EXE" == "" ]
then
	# if we could not find the MEX path, we add it manually
	MEX_EXE=/usr/local/matlab/bin/mex
fi
CFLAGS="-w -c -O3 -fopenmp" #$(OPENMP)

# Compile and link

$GCC -fPIC $CFLAGS -I$MEXPATH $INCLUDE_PATH mex_predict.cpp
$MEX_EXE mex_predict.o LDFLAGS="\$LDFLAGS -fopenmp" -lgcc ${LIBS} -outdir bin/ $MEX_ARG


