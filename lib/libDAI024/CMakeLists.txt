cmake_minimum_required(VERSION 2.4)

PROJECT(dai C CXX)

if(WIN32)
        message(STATUS "flags ${CMAKE_CXX_FLAGS}")
    SET(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /DWINDOWS /DNOMINMAX /w /EHsc /GR /MP /openmp ")
        message(STATUS "flags ${CMAKE_CXX_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_AMD64_ ")
elseif(APPLE)
    SET(CMAKE_CXX_FLAGS "-Wall -O3")
elseif(UNIX)
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -fopenmp")
    set(CMAKE_CXX_FLAGS         "${CMAKE_CXX_FLAGS} -g -pg")
    set(CMAKE_CXX_FLAGS_DEBUG   "")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
else()
    message(WARNING "unrecognized architecture")
endif()

SET(CMAKE_C_FLAGS ${CMAKE_CXX_FLAGS})

find_package(Boost 1.37.0 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

set(DAI_SRC_FILES
src/alldai.cpp
src/bp.cpp  
src/cbp.cpp 
src/emalg.cpp     
src/exceptions.cpp   
src/fbp.cpp    
src/hak.cpp    
src/properties.cpp   
src/trwbp.cpp   
src/weightedgraph.cpp
src/bbp.cpp       
src/clustergraph.cpp  
src/evidence.cpp  
src/factor.cpp       
src/gibbs.cpp  
src/jtree.cpp  
src/mf.cpp  
src/regiongraph.cpp  
src/util.cpp
src/bipgraph.cpp  
src/bp_dual.cpp  
src/daialg.cpp        
src/exactinf.cpp  
src/factorgraph.cpp  
src/graph.cpp  
src/lc.cpp     
src/mr.cpp  
src/treeep.cpp       
src/varset.cpp)        

include_directories(${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/include/dai)

message(STATUS "include: ${CMAKE_SOURCE_DIR}/include/dai")
ADD_LIBRARY(dai ${DAI_SRC_FILES})
