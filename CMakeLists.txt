cmake_minimum_required(VERSION 2.8)

project(alignement)
OPTION (SUPER4PCS_COMPILE_TESTS FALSE)
OPTION (SUPER4PCS_COMPILE_DEMOS FALSE)
add_subdirectory(ext/Super4PCS ${CMAKE_CURRENT_BINARY_DIR}/super4PCS)
add_subdirectory(projects/alignement ${CMAKE_CURRENT_BINARY_DIR}/alignement)

project(warping)
add_subdirectory(ext/ISCD_Commons ${CMAKE_CURRENT_BINARY_DIR}/commons)
add_subdirectory(ext/ISCD_LinearElasticity ${CMAKE_CURRENT_BINARY_DIR}/linearElasticity)
add_subdirectory(projects/warping ${CMAKE_CURRENT_BINARY_DIR}/warping)

project(morphing)
add_subdirectory(ext/ISCD_MshDist ${CMAKE_CURRENT_BINARY_DIR}/mshdist)
add_subdirectory(ext/tetgen ${CMAKE_CURRENT_BINARY_DIR}/tetgen)
add_subdirectory(projects/morphing ${CMAKE_CURRENT_BINARY_DIR}/morphing)

#
