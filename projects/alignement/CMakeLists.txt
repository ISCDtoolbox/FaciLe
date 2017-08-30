cmake_minimum_required(VERSION 3.0)
project(alignement)

#Super4PCS
option(SUPER4PCS_COMPILE_TESTS "" FALSE)
option(SUPER4PCS_COMPILE_DEMOS "" FALSE)
option(BUILD_SHARED_LIBS "" TRUE)
add_subdirectory(${PROJECT_SOURCE_DIR}/ext/Super4PCS ${CMAKE_CURRENT_BINARY_DIR}/super4PCS)
include_directories(${Super4PCS_INCLUDE_DIR})
link_directories(${Super4PCS_LIB_DIR})

#Eigen
find_package( Eigen3 REQUIRED )
include_directories( ${EIGEN3_INCLUDE_DIR} )

#ICP
option(BUILD_SHARED_LIBS "" TRUE)
configure_file(${PROJECT_SOURCE_DIR}/ext/CMakeLists.libICP.txt ${PROJECT_SOURCE_DIR}/ext/libICP/CMakeLists.txt)
add_subdirectory(${PROJECT_SOURCE_DIR}/ext/libICP ${CMAKE_CURRENT_BINARY_DIR}/libICP)
include_directories(${PROJECT_SOURCE_DIR}/ext/libICP/src)
link_directories(${CMAKE_BINARY_DIR}/ext/libICP)

#Commons
add_subdirectory(${PROJECT_SOURCE_DIR}/ext/Commons ${CMAKE_CURRENT_BINARY_DIR}/Commons)
include_directories(${PROJECT_SOURCE_DIR}/ext/Commons/sources)
link_directories(${CMAKE_BINARY_DIR}/ext/Commons)

# Custom definitions
add_definitions(-std=c++11 -g)

# Files
include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/src/lib)
file(GLOB LIBRARY_HEADERS ${PROJECT_SOURCE_DIR}/src/lib/*.h)
file(GLOB LIBRARY_SOURCES ${PROJECT_SOURCE_DIR}/src/lib/*.cpp)

# Library
add_library(libAlign SHARED ${LIBRARY_HEADERS} ${LIBRARY_SOURCES})
target_link_libraries(libAlign Commons super4pcs_algo super4pcs_utils libICP)
set_target_properties(libAlign PROPERTIES PREFIX "")

# Executable
add_executable(align ${PROJECT_SOURCE_DIR}/src/main.cpp)
target_link_libraries(align libAlign)
