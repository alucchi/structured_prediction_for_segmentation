cmake_minimum_required(VERSION 2.4)

PROJECT(sliceMe C CXX)

set(SLICEME_DIR ${CMAKE_SOURCE_DIR})

include(CMakeLists_common.txt)

####################################################################### LIBRARY

add_library(sliceme
${SLICEME_FILES}
)

TARGET_LINK_LIBRARIES(sliceme ${OPENCV_LIBRARIES} ${SUPERPIXEL_LIBRARY})
if(UNIX)
SET_TARGET_PROPERTIES(sliceme PROPERTIES COMPILE_FLAGS -fPIC)
endif(UNIX)
SET_TARGET_PROPERTIES(sliceme PROPERTIES ENABLE_EXPORTS TRUE)


###################################################################### BINARIES

SET_SOURCE_FILES_PROPERTIES(${SLICEME_DIR}/core/train.c PROPERTIES LANGUAGE CXX )

ADD_EXECUTABLE(train
${SLICEME_DIR}/core/train.c
${INFERENCE_FILES}
${SLICEME_FILES}
)
TARGET_LINK_LIBRARIES(train ${SLICEME_THIRD_PARTY_LIBRARIES})

ADD_EXECUTABLE(predict
${SLICEME_DIR}/core/predict.cpp
${INFERENCE_FILES}
${SLICEME_FILES}
)
TARGET_LINK_LIBRARIES(predict ${SLICEME_THIRD_PARTY_LIBRARIES})

