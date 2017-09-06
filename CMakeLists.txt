cmake_minimum_required(VERSION 2.8)
include(ExternalProject)

#Alignement project
add_subdirectory(projects/alignement ${CMAKE_CURRENT_BINARY_DIR}/alignement)

#Commons library
#add_subdirectory(${PROJECT_SOURCE_DIR}/ext/ISCD_Commons ${CMAKE_CURRENT_BINARY_DIR}/commons)
add_subdirectory(projects/morphing ${CMAKE_CURRENT_BINARY_DIR}/morphing)
add_subdirectory(projects/warping  ${CMAKE_CURRENT_BINARY_DIR}/warping)

#project(warping)
#add_subdirectory(ext/ISCD_LinearElasticity ${CMAKE_CURRENT_BINARY_DIR}/linearElasticity)
#add_subdirectory(projects/warping ${CMAKE_CURRENT_BINARY_DIR}/warping)

#project(morphing)
#add_subdirectory(ext/ISCD_MshDist ${CMAKE_CURRENT_BINARY_DIR}/mshdist)
#add_subdirectory(ext/tetgen ${CMAKE_CURRENT_BINARY_DIR}/tetgen)
#add_subdirectory(projects/morphing ${CMAKE_CURRENT_BINARY_DIR}/morphing)

find_program(TETGEN tetgen)
if(TETGEN)
  message("-- Found tetgen: ${TETGEN}")
  set(TETGEN_EXE ${TETGEN})
else()
  message("-- Will install tetgen")
  ExternalProject_Add(tetgen GIT_REPOSITORY https://github.com/ufz/tetgen.git INSTALL_COMMAND "" PREFIX ext/tetgen)
  set(TETGEN_EXE ${CMAKE_CURRENT_BINARY_DIR}/ext/tetgen/src/tetgen-build/tetgen)
endif()

find_program(MMGS mmgs_O3)
find_program(MMG3D mmgs_O3)
if(MMGS AND MMG3D)
  message("-- Found mmgs: ${MMGS}")
  message("-- Found mmg3d: ${MMG3D}")
  set(MMGS_EXE ${MMGS})
  set(MMG3D_EXE ${MMG3D})
else()
  message("-- Will install mmgTools")
  ExternalProject_Add(mmgTools GIT_REPOSITORY https://github.com/MmgTools/mmg.git INSTALL_COMMAND "" PREFIX ext/mmgTools)
  set(MMGS_EXE ${CMAKE_CURRENT_BINARY_DIR}/ext/mmgTools/src/mmgTools/bin/mmgs_O3)
  set(MMG3D_EXE ${CMAKE_CURRENT_BINARY_DIR}/ext/mmgTools/src/mmgTools/bin/mmg3d_O3)
endif()

find_program(MSHDIST mshdist)
if(MSHDIST)
  message("-- Found mshdist: ${MSHDIST}")
  set(MSHDIST_EXE ${MSHDIST})
else()
  message("-- Will install mshdist")
  ExternalProject_Add(mshdist GIT_REPOSITORY https://github.com/ISCDtoolbox/Mshdist.git INSTALL_COMMAND "" PREFIX ext/mshdist)
  set(MSHDIST_EXE ${CMAKE_CURRENT_BINARY_DIR}/ext/mshdist/src/mshdist-build/mshdist)
endif()

#Set the paths to the different executables
#Exterior tools
set(SUPER4PCS_EXE   ${CMAKE_CURRENT_BINARY_DIR}/super4PCS/demos/Super4PCS/Super4PCS )
#Compiled applications
set(MORPHING_EXE    ${CMAKE_CURRENT_BINARY_DIR}/morphing/morphing)
set(WARPING_EXE     ${CMAKE_CURRENT_BINARY_DIR}/warping/warping)
set(ALIGN_EXE       ${CMAKE_CURRENT_BINARY_DIR}/alignement/align)
#Python scripts
set(BOUNDING_EXE    "python ${CMAKE_SOURCE_DIR}/projects/tools/fastBoundingMesh.py")
set(SHELL_EXE       "python ${CMAKE_SOURCE_DIR}/projects/tools/createShell.py")
set(PYTHONICP_EXE   "python ${CMAKE_SOURCE_DIR}/projects/alignement/icp.py")
set(APPLYMATRIX_EXE "python ${CMAKE_SOURCE_DIR}/projects/tools/applyTransformationMatrix.py")
set(PCA_EXE         "python ${CMAKE_SOURCE_DIR}/projects/ACP/PCA.py")
set(FILLTOMORPH_EXE "python ${CMAKE_SOURCE_DIR}/projects/tools/fillForMorphing.py")

#Creating the python file for paths
configure_file(executable_paths.py.in ${CMAKE_SOURCE_DIR}/pipeline/executable_paths.py)
