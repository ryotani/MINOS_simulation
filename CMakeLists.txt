# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
#
#cmake -DGeant4_DIR=/home/you/geant4-install/lib64/Geant4-10.0.0 $HOME/B1
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(MINOS_sim)

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_CURRENT_SOURCE_DIR}")
find_package(ROOT 5.20 REQUIRED QUIET)
include_directories(${ROOT_INCLUDE_DIR})
#link_directories(/Users/acorsi/codes/MINOS_simulation/lib)
link_directories(/home/local1/workspace/MINOS_simulation/lib)
#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(MINOS_sim MINOS_sim.cc ${sources} ${headers})
target_link_libraries(MINOS_sim ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} ExN03BeamIn ExN03Datai ExN03Setup)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build betaPlus. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#


