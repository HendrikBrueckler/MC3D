### Find dependencies
# gmp
find_package(GMP REQUIRED)
find_package(GMPXX REQUIRED)
list(APPEND MC3D_LIB_LIST ${GMP_LIBRARIES} ${GMPXX_LIBRARY})

# OVM
if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    set(CMAKE_CXX_STANDARD 17)
    set(CXX_STANDARD_REQUIRED ON)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh/CMakeLists.txt")
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh extern/OpenVolumeMesh EXCLUDE_FROM_ALL)
    else()
        find_package(OpenVolumeMesh REQUIRED)
    endif()
endif()
list(APPEND MC3D_LIB_LIST OpenVolumeMesh::OpenVolumeMesh)

# glog
if (NOT TARGET glog::glog)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/glog/CMakeLists.txt")
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/glog extern/glog EXCLUDE_FROM_ALL)
    else()
        find_package(glog REQUIRED)
    endif()
endif()
list(APPEND MC3D_LIB_LIST glog::glog)
if (NOT MC3D_ENABLE_LOGGING)
    list(APPEND MC3D_COMPILE_DEFINITIONS_PRV "GOOGLE_STRIP_LOG=10")
endif()

# TrulySeamless3D
if (NOT TARGET TS3D::TS3D)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/TrulySeamless3D/CMakeLists.txt")
        if (CMAKE_BUILD_TYPE EQUAL "DEBUG")
            option(TRULYSEAMLESS_ENABLE_LOGGING "" ON)
        else()
            option(TRULYSEAMLESS_ENABLE_LOGGING "" OFF)
        endif()
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/TrulySeamless3D extern/TrulySeamless3D EXCLUDE_FROM_ALL)
    else()
        find_package(TS3D REQUIRED)
    endif()
endif()
list(APPEND MC3D_LIB_LIST TS3D::TS3D)

### Source files
list(APPEND MC3D_SOURCE_LIST
     "Algorithm/SingularityInitializer.cpp"
     "Algorithm/MCBuilder.cpp"
     "Algorithm/MCReducer.cpp"
     "Algorithm/MotorcycleSpawner.cpp"
     "Algorithm/MotorcycleTracer.cpp"
     "Data/Motorcycle.cpp"
     "Data/Transition.cpp"
     "Data/UVWDir.cpp"
     "Data/NodeCoordination.cpp"
     "Interface/MCGenerator.cpp"
     "Interface/Reader.cpp"
     "Interface/Writer.cpp"
     "Mesh/MCMeshManipulator.cpp"
     "Mesh/MCMeshNavigator.cpp"
     "Mesh/MCMeshProps.cpp"
     "Mesh/TetMeshManipulator.cpp"
     "Mesh/TetMeshNavigator.cpp"
     "Mesh/TetMeshProps.cpp")

### Create target
add_library(MC3D STATIC ${MC3D_SOURCE_LIST})
add_library(MC3D::MC3D ALIAS MC3D)

### Configure install paths
include(CMakePackageConfigHelpers)
if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR "CMake")
else()
  set(DEF_INSTALL_CMAKE_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/MC3D")
endif()
set(MC3D_INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH "Installation directory for CMake files")
set(MC3D_INSTALL_LIB_DIR "${CMAKE_INSTALL_LIBDIR}" CACHE PATH "Installation directory for libraries")
set(MC3D_INSTALL_INCLUDE_DIR "${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Installation directory for header files")
# ... for the build tree
set(CONF_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/../include")
set(CONF_LIB_DIRS "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
configure_package_config_file(../MC3DConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/MC3DConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")
# ... for the install tree
set(CONF_INCLUDE_DIRS "${MC3D_INSTALL_INCLUDE_DIR}")
set(CONF_LIB_DIRS "${MC3D_INSTALL_LIB_DIR}")
configure_package_config_file(../MC3DConfig.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/MC3DConfig.cmake"
  INSTALL_DESTINATION "${MC3D_INSTALL_CMAKE_DIR}")

### Include own headers with public access
target_include_directories(MC3D
                           PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../include>"
                           "$<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>")

### Link with dependencies
target_link_libraries_system(MC3D PUBLIC ${MC3D_LIB_LIST})

### Set compile options
# properties
set_target_properties(MC3D PROPERTIES
                      DEBUG_POSTFIX "d"
                      VERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
                      SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}
                      CXX_STANDARD 17
                      CXX_STANDARD_REQUIRED ON)
# compile options
list(APPEND MC3D_COMPILE_OPTIONS_PRV "-Wall" "-Wextra" "-Werror" "-Wshadow" "-Wunused" "-Wunused-function")
target_compile_options(MC3D PRIVATE ${MC3D_COMPILE_OPTIONS_PRV})# preprocessor defines
# preprocessor defines
target_compile_definitions(MC3D PRIVATE ${MC3D_COMPILE_DEFINITIONS_PRV})

### Include headers as system headers to prevent warnings/errors from strict compile definitions
# target_include_directories(MC3D SYSTEM PUBLIC ${GMPXX_INCLUDE_DIR})
# target_include_directories(MC3D SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh/src)
# target_include_directories(MC3D SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/glog/src)
# target_include_directories(MC3D SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/TrulySeamles3D/include)

### Install and export
install(DIRECTORY "${PROJECT_SOURCE_DIR}/include/MC3D"
        DESTINATION "${MC3D_INSTALL_INCLUDE_DIR}")

install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/MC3DConfig.cmake"
        DESTINATION "${MC3D_INSTALL_CMAKE_DIR}")

install(TARGETS       MC3D
        EXPORT        MC3DTargets
        LIBRARY       DESTINATION "${MC3D_INSTALL_LIB_DIR}"
        ARCHIVE       DESTINATION "${MC3D_INSTALL_LIB_DIR}")

install(EXPORT MC3DTargets
        FILE MC3DTargets.cmake
        NAMESPACE MC3D::
        DESTINATION "${MC3D_INSTALL_CMAKE_DIR}")

export(EXPORT MC3DTargets
       NAMESPACE MC3D::)
