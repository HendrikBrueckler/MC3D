find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
    # Update submodules as needed
    message(STATUS "Submodule init")
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive -f
                    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                    RESULT_VARIABLE GIT_SUBMOD_RESULT)
    if(NOT GIT_SUBMOD_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
    endif()
    message(STATUS "Submodule patch")
    execute_process(COMMAND ${GIT_EXECUTABLE} apply "../patches/OpenVolumeMesh.patch"
                    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/"
                    RESULT_VARIABLE GIT_PATCH_RESULT)
    if(NOT GIT_PATCH_RESULT EQUAL "0")
        message(FATAL_ERROR "git apply ../patches/OpenVolumeMesh.patch' failed with ${GIT_PATCH_RESULT}")
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh/CMakeLists.txt"
   OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/googletest/CMakeLists.txt"
   OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/glog/CMakeLists.txt"
   OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/CLI11/CMakeLists.txt")
    message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()
