
---------------------------------------------------------------------------------

*** Acknowledgment ***

This code was developped by Aurelien Lucchi during his doctoral studies in the computer vision laboratory run by Prof. Pascal Fua at EPFL.
http://cvlab.epfl.ch/

Contributors:
- Aurelien Lucchi
- Yunpeng Li
- Pol Monsó Purtí

---------------------------------------------------------------------------------

*** How to compile the code ***

This tutorial will help you compile the code and run it on your own system. It is intended for people with a Linux system but windows or mac users should also be able to compile the code.

Step-by-step instructions:

1. Copy the source code from the web.

2. Make sure you have cmake installed on your system (cmake is a tool to create makefiles).

3. Compile third-party libraries :
[UNIX]
- OpenCv (Install package opencv-dev in Ubuntu or download from http://opencv.org)
- zlib (http://www.zlib.net/)
- SLIC: Go to lib/slic and type "cmake .; make"
- ITK : Download from the web site and set review flag to ON (can be done by using advanced mode with ccmake and settting variable USE_REVIEW to ON).
- Go to lib/libDAI-0.2.4 and type "make". Note that libDAI requires the following dependencies for Ubuntu:
 sudo apt-get install g++ make doxygen graphviz libboost-dev libboost-graph-dev libboost-program-options-dev

[WINDOWS]
- OpenCV: download archive from http://opencv.org. Then create makefile with cmake and compile with "MSBuild.exe opencv.sln"
- SLIC: Go to lib/slic and type "cmake .; make"
- ITK : Download from the web site and set review flag to ON (can be done by using advanced mode with ccmake and settting variable USE_REVIEW to ON).
- Go to lib/libDAI-0.2.4 and type "cmake" and "MSBuild.exe dai.sln". Note that libDAI requires boost.
- [OPTIONAL] Install 7-zip to compress intermediary results

4. Main code
You can edit CMakeLists_common.txt or just change the flags with the ccmake interface.
You might want to turn off some of the dependencies. Look at the USE_??? flags.

Go to the superpixels_public directory and compile with
[UNIX]
cmake .; make
[WINDOWS]
"cmake ." and "MSBuild.exe sliceMe.sln" -or- open the sliceMe.sln with Visual Studio.

5. Edit the config file

This program relies on a configuration file that contains all the parameters needed for training or applying the algorithm to a new dataset. Below is a description of each parameter.

trainingDir: Input volume (TIF cube).
maskTrainingDir: Mask volume (TIF cube or series of png images containing black and white voxels only)
outputModel: Output model file (default=model.txt)
testDir: Test input volume (TIF cube).
maskTestDir: Test mask volume (TIF cube containing black and white voxels only).
stepForOutputFiles: output files and corresponding scores are computed every stepForOutputFiles iterations. Reasonable values are 10 for default training algorithm and 100 for stochastic optimization (-w 9 option described below).
featureTypes: Type of features extracted from the image/cube and used to train a model. The default value is featureTypes=33.
nMaxIterations: Maximum number of iterations. This number should be around 500-1000 iterations depending on how much training data is provided to the algorithm.

An example configuration file named config.txt with a set of default parameters is provided in the root directory or on the web at http://cvlab.epfl.ch/software/biomedplugins.

6. Running the code
Look at howto.pdf

---------------------------------------------------------------------------------

*** Additional information for developers ***

- Multithreading
For ssvm, edit svm_struct_globals.h and change NTHREADS

