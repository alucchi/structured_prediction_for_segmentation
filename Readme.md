
---------------------------------------------------------------------------------

# Acknowledgment #

This code was developped by Aurelien Lucchi during his doctoral studies in the [computer vision laboratory](http://cvlab.epfl.ch/) run by Prof. Pascal Fua at EPFL.


Contributors:

- Aurelien Lucchi (github [@alucchi](https://github.com/alucchi))
- Yunpeng Li
- Pol Monsó Purtí (github [@quimnuss](https://github.com/quimnuss))

---------------------------------------------------------------------------------

# How to compile the code #

This tutorial will help you compile the code and run it on your own system. It is intended for people with a Linux system but windows or mac users should also be able to compile the code.

Step-by-step instructions:

1. Copy the source code from the web.

2. Make sure you have cmake installed on your system (cmake is a tool to create makefiles).

3. Compile third-party libraries :

	**UNIX**
	
	- OpenCv (Install package opencv-dev in Ubuntu or download from [http://opencv.org](http://opencv.org)). Versions 3 and 2.4 work.
	- [zlib](http://www.zlib.net/)
	
	If you have problems with zlib, we recommend using the one shipped as a 3rdparty within `opencv`. Set `ZLIB_INCLUDE_DIR` to `3rdparty/zlib` and `ZLIB_LIBRARY` to `build/3rdpartz/lib/zlib.lib`
	
	- SLIC: Go to `lib/slic` and type `cmake .; make`
	- ITK : Download from the web site and set review flag to `ON` (can be done by using advanced mode with ccmake and settting variable `Module_ITKReview` to ON). Tested version 4.7.
	- Go to `lib/libDAI-0.2.4` and type `make`. Note that `libDAI` requires the following dependencies for Ubuntu:
	 `sudo apt-get install g++ make doxygen graphviz libboost-dev libboost-graph-dev libboost-program-options-dev`
	 
	**WINDOWS**

    We recommend using `MSBuild` to build. For example `"C:\Windows\Microsoft.NET\Framework64\v4.0.30319\MSBuild.exe" opencv.sln /p:Configuration=Release /m`
	
	- [OpenCV](http://opencv.org): download archive from web site. Then create makefile with cmake and compile with `MSBuild.exe opencv.sln`. Versions 3 and 2.4 work.
	- [zlib](http://www.zlib.net/)
	
	If you have problems with zlib, we recommend using the one shipped as a 3rdparty within `opencv`. Set `ZLIB_INCLUDE_DIR` to `3rdparty/zlib` and `ZLIB_LIBRARY` to `build/3rdparty/lib/Release/zlib.lib`. In OpenCV 2.4 you may have to copy the generated `zconf.h` from `build/3rdparty/include` or `sources/build/3rdparty/zlib` to `3rdparty/zlib`, or include both paths.
	
	- SLIC: Go to `lib/slic` and type `cmake .; make`
	- [ITK](http://itk.org): Download from the web site and set review flag to `ON` (can be done by using advanced mode with ccmake and settting variable `Module_ITKReview` to `ON`). Tested version 4.7.
	- Go to `lib/libDAI-0.2.4` and type `cmake` and `MSBuild.exe dai.sln`. Note that `libDAI` requires `boost`.
    - Add [`dirent.h`](http://www.softagalleria.net/download/dirent/) to `core/`
	- [OPTIONAL] Install `7-zip` to compress intermediary results.


4. Main code

	You can edit `CMakeLists_common.txt` or just change the flags with the `ccmake` (UNIX) or `cmake-gui` (WINDOWS) interface.
	You might want to turn off some of the dependencies. Look at the `USE_???` flags.

	Go to the `superpixels_public` directory and compile with
	
	**UNIX**
	
	`mkdir build
	ccmake ..`
	
	Once configured and generated run
	`make`
	
	**WINDOWS**
		
	`mkdir build
	cmake-gui ..`

	Once configured and generated run
	`MSBuild.exe sliceMe.sln` _or_ open the `sliceMe.sln` with Visual Studio.

5. Edit the config file

	This program relies on a configuration file that contains all the parameters needed for training or applying the algorithm to a new dataset. Below is a description of each parameter.

	* `trainingDir`: Input volume (TIF cube).
	* `maskTrainingDir`: Mask volume (TIF cube or series of png images containing black and white voxels only)
	* `outputModel`: Output model file (default=`model.txt`)
	* `testDir`: Test input volume (TIF cube).
	* `maskTestDir`: Test mask volume (TIF cube containing black and white voxels only).
	* `stepForOutputFiles`: output files and corresponding scores are computed every stepForOutputFiles iterations. Reasonable values are 10 for default training algorithm and 100 for stochastic optimization (`-w 9` option described below).
	* `featureTypes`: Type of features extracted from the image/cube and used to train a model. The default value is `featureTypes=33`.
	* `nMaxIterations`: Maximum number of iterations. This number should be around 500-1000 iterations depending on how much training data is provided to the algorithm.

	An example configuration file named `config.txt` with a set of default parameters is provided in the root directory or on the [biomed web](http://cvlab.epfl.ch/software/biomedplugins).

6. Running the code

	Look at `howto.pdf`

---------------------------------------------------------------------------------

## Additional information for developers ##

- Multithreading

For ssvm, edit `svm_struct_globals.h` and change `NTHREADS`

