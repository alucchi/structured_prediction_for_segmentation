
*******************************************************************************
How to compile the mex files?

Your first have to compile the mex code using the file in matlab/mex/buildmex.m. You have to edit the file and change the path to point to the libraries on your system.
Once you got the code to compile, you can test it with the file in matlab/src/test_predict.m

For Linux users:
You can also use the shell script mex/buildmex.sh. You need to specify the path of the developement libraries used to compile the code. You can then run the script by typing:
 bash buildmex.sh

*******************************************************************************
How to start matlab?
LD_LIBRARY_PATH=.:$PATH:$LD_LIBRARY_PATH LD_PRELOAD=/usr/lib/gcc/x86_64-linux-gnu/4.6/libstdc++.so /usr/local/MATLAB/R2013a/bin/matlab

*******************************************************************************
