cmake_minimum_required(VERSION 2.4)

PROJECT(sliceMe)

set(SLICEME_DIR ${CMAKE_SOURCE_DIR}/../)
include(${SLICEME_DIR}CMakeLists_common.txt)

####################################################################### LIBRARY


###################################################################### BINARIES

ADD_EXECUTABLE(exportFeatures
exportFeatures.cpp
${INFERENCE_FILES}
${SLICEME_FILES}
)
TARGET_LINK_LIBRARIES(exportFeatures ${SLICEME_THIRD_PARTY_LIBRARIES})

