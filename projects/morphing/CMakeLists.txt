cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
set(CMAKE_C_FLAGS   "-w -g -O3")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#project(libMorphing)
file(                  GLOB_RECURSE sources src/*)
file(                  GLOB_RECURSE headers src/*.h)
add_library(           libMorphing SHARED ${sources})
set_target_properties( libMorphing PROPERTIES PREFIX "")

#target_link_libraries( libMorphing Commons)
#
#find_library(          Elasticity NAMES Elasticity HINTS "/usr/local/lib/iscd" required)
#link_directories(      "/usr/local/lib/iscd")
#include_directories(   "/usr/local/include/iscd")
if(${Commons})
  message("toto")
  target_link_libraries( libMorphing Commons)
else()
  find_library(          Commons NAMES Commons HINTS $ENV{HOME} required)
  target_link_libraries( libMorphing Commons)
endif()



#project(morphing)
#file(                  GLOB_RECURSE sources src/main.c)
add_executable(        morphing src/main.c)
target_link_libraries( morphing libMorphing ${Commons})

#install(               TARGETS libMorphing LIBRARY DESTINATION "/usr/local/lib/")
#install(               FILES   ${headers} DESTINATION "/usr/local/include/morphing/")
#set_property(          TARGET  libMorphing PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
#set_property(          TARGET  libMorphing PROPERTY OUTPUT_NAME Morphing)
#install(               TARGETS morphing RUNTIME  DESTINATION "/usr/local/bin/")
#set_property( TARGET morphing PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
