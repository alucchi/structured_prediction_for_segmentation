cmake_minimum_required(VERSION 2.4)

PROJECT(supervoxel)
INCLUDE_REGULAR_EXPRESSION("^.*$")


#SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE Release)

if(WIN32)
    SET(CMAKE_CXX_FLAGS "/DWIN32 /D_WINDOWS /W3 /Zi /EHsc /GR /MP /openmp")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /D_AMD64_")
elseif(APPLE)
    SET(CMAKE_CXX_FLAGS "-Wall -O3")
elseif(UNIX)
    set(CMAKE_CXX_FLAGS         "-g -pg")
    set(CMAKE_CXX_FLAGS_DEBUG   "")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
else()
    message(WARNING "unrecognized architecture")
endif()

set(CMAKE_C_FLAGS         ${CMAKE_CXX_FLAGS})
set(CMAKE_C_FLAGS_DEBUG   ${CMAKE_CXX_FLAGS_DEBUG})
set(CMAKE_C_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})

set(OpenCV_DIR ~/src/third_party/opencv-2.4.9/)
#set(OpenCV_DIR ~/src/ThirdParty/OpenCV-2.4.2/)
FIND_PACKAGE(OpenCV REQUIRED)
# add include directories
INCLUDE_DIRECTORIES(${OPENCV_INCLUDE_DIR})

add_library(supervoxel
LKM.cpp
utils.cpp)
TARGET_LINK_LIBRARIES(supervoxel ${OpenCV_LIBS})
SET_TARGET_PROPERTIES(supervoxel PROPERTIES COMPILE_FLAGS -fPIC)

#ADD_EXECUTABLE(supervoxel_test
#supervoxel_test.cpp
#LKM.cpp
#utils.cpp)
#TARGET_LINK_LIBRARIES(supervoxel_test ${OpenCV_LIBS})
#
#ADD_EXECUTABLE(superpixel_test
#superpixel_test.cpp
#LKM.cpp
#utils.cpp)
#TARGET_LINK_LIBRARIES(superpixel_test ${OpenCV_LIBS})
