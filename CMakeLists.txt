project(USRPutils C CXX)
cmake_minimum_required(VERSION 2.8)

set(USRPutils_MAJOR_VERSION 0)
set(USRPutils_MINOR_VERSION 3)
set(USRPutils_PATCH_VERSION 0)
set(USRPutils_VERSION ${USRPutils_MAJOR_VERSION}.${USRPutils_MINOR_VERSION}.${USRPutils_PATH_VERSION})


#
#Set General Compile Flags
#
set(CMAKE_CXX_FLAGS "-O2 -march=native")
set(CMAKE_C_FLAGS "-O2 -march=native")

if(DEBUG)
  set(CMAKE_CXX_FLAGS "-DDEBUG -g ${CMAKE_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "-DDEBUG -g ${CMAKE_CXX_FLAGS}")
endif()

#Set the runtime executable output:
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${USRPutils_SOURCE_DIR}/build)
#
#
#Check for required packages
#
set(CMAKE_MODULE_PATH ${USRPutils_SOURCE_DIR}/cmake/Modules)
find_package(UHD REQUIRED)
find_package(FFTW3f REQUIRED)
find_package(Threads REQUIRED)
find_package(Boost COMPONENTS system REQUIRED)

add_subdirectory(src)
