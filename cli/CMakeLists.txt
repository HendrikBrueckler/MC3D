set(CLI_NAME "mc3d_cli")

if (NOT TARGET CLI11::CLI11)
    if(EXISTS "${PROJECT_SOURCE_DIR}/extern/CLI11/CMakeLists.txt")
        add_subdirectory(${PROJECT_SOURCE_DIR}/extern/CLI11 extern/CLI11 EXCLUDE_FROM_ALL)
    else()
        find_package(CLI11 REQUIRED)
    endif()
endif()

add_executable(${CLI_NAME} main.cpp)
target_link_libraries(${CLI_NAME} MC3D::MC3D CLI11::CLI11 OpenVolumeMesh::OpenVolumeMesh)
