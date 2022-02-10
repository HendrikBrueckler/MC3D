option(MC3D_SUBMODULES_MANUAL "Download MC3D submodules manually" OFF)

find_package(Git QUIET)

set(OVM_PATCH_TIMESTAMP "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/ovmpatch.stamp")
if(NOT ${MC3D_SUBMODULES_MANUAL})
    if(GIT_FOUND AND EXISTS "${CMAKE_SOURCE_DIR}/.git")
        git_update_submodules("${CMAKE_CURRENT_SOURCE_DIR}")
        git_patch("../patches/OpenVolumeMesh.patch" "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/")
        execute_process(COMMAND ${CMAKE_COMMAND} -E touch ${OVM_PATCH_TIMESTAMP})
    else()
        message(WARNING "TRULYSEAMLESS_SUBMODULES_MANUAL is on, but this is not a git repository")
    endif()
else()
    if (GIT_FOUND AND NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/ovmpatch.stamp")
        git_patch("../patches/OpenVolumeMesh.patch" "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/")
        execute_process(COMMAND ${CMAKE_COMMAND} -E touch ${OVM_PATCH_TIMESTAMP})
    endif()
    if (NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/ovmpatch.stamp")
        message(WARNING "Could not patch OVM, proceeding with unpatched version, AVOID SELFADJACENT BLOCKS")
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/OpenVolumeMesh/CMakeLists.txt"
   OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/googletest/CMakeLists.txt"
   OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/glog/CMakeLists.txt"
   OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/CLI11/CMakeLists.txt"
   OR NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/TrulySeamless3D/CMakeLists.txt")
    if (${MC3D_STANDALONE})
        message(FATAL_ERROR "The submodules were not downloaded! Please update submodules and try again.")
    else()
        message(WARNING "MC3D submodules not downloaded! Must supply them in a base project.")
    endif()
endif()
