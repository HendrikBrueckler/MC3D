include(MC3DMacros)

### Better safe than sorry...
list(APPEND MC3D_LIB_COMPILE_OPTIONS "-Wall" "-Wextra" "-Werror" "-Wshadow")

### Add all source files
list(APPEND SOURCE_LIST
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
     "Mesh/TetMeshProps.cpp"
     )
add_library(${MC3D_LIB_NAME} STATIC ${SOURCE_LIST})

### Find dependencies
# gmp
find_package(GMP REQUIRED)
find_package(GMPXX REQUIRED)
list(APPEND LIB_LIST ${GMP_LIBRARIES} ${GMPXX_LIBRARY})

# OVM
add_subdirectory(${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh extern/OpenVolumeMesh EXCLUDE_FROM_ALL)
list(APPEND LIB_LIST OpenVolumeMesh::OpenVolumeMesh)

# glog
add_subdirectory(${PROJECT_SOURCE_DIR}/extern/glog extern/glog EXCLUDE_FROM_ALL)
list(APPEND LIB_LIST glog::glog)
if (NOT MC3D_ENABLE_LOGGING)
    list(APPEND MC3D_LIB_COMPILE_DEFINITIONS "GOOGLE_STRIP_LOG=10")
endif()

### Link with dependencies
target_link_libraries_system(${MC3D_LIB_NAME} PUBLIC ${LIB_LIST})

### Set compile options
set_target_properties(${MC3D_LIB_NAME} PROPERTIES CXX_STANDARD 11 CXX_STANDARD_REQUIRED ON)
target_compile_options(${MC3D_LIB_NAME} PRIVATE ${MC3D_LIB_COMPILE_OPTIONS})
target_compile_definitions(${MC3D_LIB_NAME} PRIVATE ${MC3D_LIB_COMPILE_DEFINITIONS})

### Include own headers with public access
target_include_directories(${MC3D_LIB_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

### Include OVM, gmp, glog headers as system headers to prevent warnings/errors from strict compile definitions
target_include_directories(${MC3D_LIB_NAME} SYSTEM PUBLIC ${GMPXX_INCLUDE_DIR})
target_include_directories(${MC3D_LIB_NAME} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh/src)
target_include_directories(${MC3D_LIB_NAME} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/extern/glog/src)
