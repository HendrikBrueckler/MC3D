cmake_minimum_required(VERSION 3.10)

project(MC3D VERSION 1.0)

### Default to C++11, CLI and tests need C++17, however
set(CMAKE_CXX_STANDARD 11)

### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build this in a source directory (or any directory with a CMakeLists.txt file). "
                        "Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

### Cmake package finder files
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### Build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE "Release" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
  "Debug" "Release" "RelWithDebInfo")
endif()

if(${CMAKE_PROJECT_NAME} STREQUAL ${PROJECT_NAME})
    set(MC3D_STANDALONE TRUE)
else()
    set(MC3D_STANDALONE FALSE)
endif()

### Binary output
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/lib/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/Build/bin/)

### Options
option(MC3D_ENABLE_LOGGING  "Enable logging for MC3D" ${MC3D_STANDALONE})
option(MC3D_BUILD_CLI       "Build CLI app for MC3D"  ${MC3D_STANDALONE})
option(MC3D_BUILD_TESTS     "Build tests for MC3D"    ${MC3D_STANDALONE})

### Load external libs
add_subdirectory(extern)

### Add library
set(MC3D_LIB_NAME mc3d)
add_subdirectory(src)

### Add CLI
if (MC3D_BUILD_CLI)
    add_subdirectory(cli)
endif()

### Tests
if (MC3D_BUILD_TESTS)
  add_subdirectory(tests)
endif()
