set(RENDERER_NAME "mc3d_renderer")

add_executable(${RENDERER_NAME} main.cpp)
target_link_libraries(${RENDERER_NAME} MC3D::MC3D)
target_include_directories(${RENDERER_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/extern/volumeshOS_core/src)
set_target_properties(${RENDERER_NAME} PROPERTIES
                        CXX_STANDARD 17
                        CXX_STANDARD_REQUIRED ON
                        RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/renderer)
