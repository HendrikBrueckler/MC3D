set(CLI_NAME "mc3d_cli")

add_subdirectory(${PROJECT_SOURCE_DIR}/extern/CLI11 extern/CLI11 EXCLUDE_FROM_ALL)

add_executable(${CLI_NAME} main.cpp)
target_link_libraries(${CLI_NAME} ${MC3D_LIB_NAME} CLI11::CLI11)
set_target_properties(${CLI_NAME} PROPERTIES
                        CXX_STANDARD 17
                        CXX_STANDARD_REQUIRED ON
                        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/cli)
