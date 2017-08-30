cmake_minimum_required(VERSION 2.8)


#Alignement project
add_subdirectory(projects/alignement ${CMAKE_CURRENT_BINARY_DIR}/alignement)

#Commons library
#add_subdirectory(${PROJECT_SOURCE_DIR}/ext/ISCD_Commons ${CMAKE_CURRENT_BINARY_DIR}/commons)


#project(warping)
#add_subdirectory(ext/ISCD_LinearElasticity ${CMAKE_CURRENT_BINARY_DIR}/linearElasticity)
#add_subdirectory(projects/warping ${CMAKE_CURRENT_BINARY_DIR}/warping)

#project(morphing)
#add_subdirectory(ext/ISCD_MshDist ${CMAKE_CURRENT_BINARY_DIR}/mshdist)
#add_subdirectory(ext/tetgen ${CMAKE_CURRENT_BINARY_DIR}/tetgen)
#add_subdirectory(projects/morphing ${CMAKE_CURRENT_BINARY_DIR}/morphing)

#Set the paths to the different executables
#Exterior tools
set(MSHDIST_EXE   ${CMAKE_CURRENT_BINARY_DIR}/mshdist/mshdist)
set(TETGEN_EXE    ${CMAKE_CURRENT_BINARY_DIR}/tetgen/tetgen)
set(MMGS_EXE      ${CMAKE_CURRENT_BINARY_DIR}/mmgTools/bin/mmgs_O3)
set(MMG3D_EXE     ${CMAKE_CURRENT_BINARY_DIR}/mmgTools/bin/mmg3d_O3)
set(SUPER4PCS_EXE ${CMAKE_CURRENT_BINARY_DIR}/super4PCS/demos/Super4PCS/Super4PCS )
#Compiled applications
set(MORPHING_EXE  ${CMAKE_CURRENT_BINARY_DIR}/morphing/morphing)
set(WARPING_EXE   ${CMAKE_CURRENT_BINARY_DIR}/warping/warping)
set(ALIGN_EXE     ${CMAKE_CURRENT_BINARY_DIR}/alignement/align)
#Python scripts
set(BOUNDING_EXE  ${CMAKE_SOURCE_DIR}/projects/tools/fastBoundingMesh.py)
#Creating the python file
configure_file(executable_paths.py.in ${CMAKE_SOURCE_DIR}/projects/executable_paths.py)
