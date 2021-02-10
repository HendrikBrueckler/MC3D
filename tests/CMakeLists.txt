include(GoogleTest)

option(INSTALL_GTEST OFF)
add_subdirectory(${PROJECT_SOURCE_DIR}/extern/googletest extern/googletest)

set(MC3D_TESTUTILS_LIB "mc3dtestutils")
add_library(${MC3D_TESTUTILS_LIB} STATIC TestUtils.cpp)
target_link_libraries(${MC3D_TESTUTILS_LIB} PRIVATE PUBLIC gtest gmock ${MC3D_LIB_NAME})
target_compile_definitions(${MC3D_TESTUTILS_LIB} PRIVATE "-DTEST_RESOURCES_PATH=\"${PROJECT_SOURCE_DIR}/tests/resources/\"")
set_target_properties(${MC3D_TESTUTILS_LIB} PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

include(MC3DTestMacros)

mc3d_add_test(BlackBoxTest BlackBoxTest.cpp)
mc3d_add_test(ReaderTest ReaderTest.cpp)
mc3d_add_test(InitializerTest InitializerTest.cpp)
mc3d_add_test(MotorcycleSpawnerTest MotorcycleSpawnerTest.cpp)
mc3d_add_test(MotorcycleTracerTest MotorcycleTracerTest.cpp)
mc3d_add_test(MCBuilderTest MCBuilderTest.cpp)
mc3d_add_test(MCReducerTest MCReducerTest.cpp)
