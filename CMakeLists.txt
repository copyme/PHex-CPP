cmake_minimum_required(VERSION 3.20)
project(PHex-CPP)

##### helpers

function(download_file url filename)
  if(NOT EXISTS ${filename})
    file(DOWNLOAD ${url} ${filename} TIMEOUT 360)
  endif()
endfunction(download_file)
#####

execute_process(COMMAND git submodule update --init -- ./extern/libigl
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

execute_process(COMMAND git submodule update --init -- ./extern/directional
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

execute_process(COMMAND git submodule update --init -- ./extern/ceres-solver
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
                   
execute_process(COMMAND git submodule update --init -- ./extern/libhedra
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

execute_process(COMMAND git submodule update --remote
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

set(BUILD_SHARED_LIBS OFF CACHE BOOL "No shared" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "No testing targets" FORCE)
   
if(UNIX AND NOT APPLE)
  set(CGAL_Boost_USE_STATIC_LIBS OFF CACHE BOOL "No static boost libs" FORCE)
endif()
   

set(LIBIGL_USE_STATIC_LIBRARY OFF CACHE BOOL "No static libs" FORCE)
set(LIBIGL_EXPORT_TARGETS OFF CACHE BOOL "No igl targets" FORCE) # exports Eigen, etc.
set(LIBIGL_WITH_EMBREE OFF CACHE BOOL "No embree" FORCE)
set(LIBIGL_WITH_COMISO OFF CACHE BOOL "No comiso" FORCE)
set(LIBIGL_WITH_PYTHON OFF CACHE BOOL "No python" FORCE)
set(LIBIGL_WITH_TETGEN OFF CACHE BOOL "No tetgen" FORCE)
set(LIBIGL_WITH_TRIANGLE OFF CACHE BOOL "No triqngle" FORCE)
set(LIBIGL_WITH_PNG OFF CACHE BOOL "No png" FORCE)
set(LIBIGL_WITH_OPENGL OFF CACHE BOOL "No OpenGL" FORCE)
set(LIBIGL_WITH_OPENGL_GLFW OFF CACHE BOOL "No GLFW" FORCE)
set(LIBIGL_WITH_OPENGL_GLFW_IMGUI OFF CACHE BOOL "No IMGUI" FORCE)
set(LIBIGL_WITH_XML OFF  CACHE BOOL "NO XML" FORCE)
set(LIBIGL_WITH_PREDICATES OFF  CACHE BOOL "No predicqtes" FORCE)
set(BUILD_EXAMPLES OFF CACHE BOOL "No examples" FORCE)
set(GFLAGS OFF CACHE BOOL "No GFLAGS" FORCE)
set(CXSPARSE OFF CACHE BOOL "No CXSPARSE" FORCE)
set(SUITESPARSE OFF CACHE BOOL "NO Eigen SuiteSparse (problems on Win)" FORCE)

set(LIBIGL_WITH_MATLAB ON CACHE BOOL "NO MATLAB" FORCE)
set(LIBIGL_WITH_CGAL ON CACHE BOOL "Use CGAL" FORCE)
set(MINIGLOG ON CACHE BOOL "Use MINIGLOG" FORCE)
set(EIGENSPARSE ON CACHE BOOL "NO Eigen sparse (ceres)" FORCE)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/libigl ${CMAKE_BINARY_DIR}/extern/libigl)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/ceres-solver ${CMAKE_BINARY_DIR}/extern/ceres-solver)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/directional/include/ ${CMAKE_BINARY_DIR}/extern/directional)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/libhedra/include/ ${CMAKE_BINARY_DIR}/extern/libhedra)


# options
option(CPF_WITH_FAST_MATH            "Compile fast_sparse from gptoobox"                   ON)

if (CPF_WITH_FAST_MATH)
  download_file(https://github.com/alecjacobson/gptoolbox/archive/refs/heads/master.zip
          ${CMAKE_CURRENT_SOURCE_DIR}/extern/gptoolbox-master.zip)
  file(ARCHIVE_EXTRACT INPUT ${CMAKE_CURRENT_SOURCE_DIR}/extern/gptoolbox-master.zip DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/extern/)
  file(REMOVE ${CMAKE_CURRENT_SOURCE_DIR}/extern/gptoolbox-master.zip)

  # set already existing igl
  set(LIBIGL_INCLUDE_DIR ${CMAKE_BINARY_DIR}/extern/libigl/include FORCE)
  set(LIBIGL_FOUND TRUE CACHE BOOL "Skip tests" FORCE)

  set(WITH_ELTOPO OFF CACHE BOOL "No ELTOP" FORCE)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/extern/gptoolbox-master/mex/)
endif()
# Workaround to be able to deal with compiler bug "Too many sections" with CYGWIN.
if( CYGWIN )
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -Wall -O3 -Wa,-mbig-obj" )
endif()

set(CPF_DIR "" CACHE PATH "Path to the CPF MATLAB files")

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Dualizer)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Planarize)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Parameterization)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Mesher)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/MeshCutter)

if (EXISTS ${CPF_DIR})
  if(CYGWIN OR MSVC)
    add_custom_target(cpf_install
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/MeshCutter/MeshCutter_bin.exe
            ${CPF_DIR}/external/MeshCutter/MeshCutter_bin.exe

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/Mesher/Mesher_bin.exe
            ${CPF_DIR}/external/Mesher/Mesher_bin.exe

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/Parametrization/Parametrization_bin.exe
            ${CPF_DIR}/external/Parametrization/Parametrization_bin.exe

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/Dualizer/Dualizer_bin.exe
            ${CPF_DIR}/external/Dualizer/Dualizer_bin.exe

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/Planarize/Planarize_bin.exe
            ${CPF_DIR}/external/Planarizer/Planarize_bin.exe

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/extern/gptoolbox-master/mex/fast_sparse.mexw64
            ${CPF_DIR}/external/fast_sparse/fast_sparse.mexw64
            )
  else()
    add_custom_target(cpf_install
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/MeshCutter/MeshCutter_bin
            ${CPF_DIR}/external/MeshCutter/MeshCutter_bin

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/Mesher/Mesher_bin
            ${CPF_DIR}/external/Mesher/Mesher_bin

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/Parameterization/Parameterization_bin
            ${CPF_DIR}/external/Parameterization/Parameterization_bin

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/Dualizer/Dualizer_bin
            ${CPF_DIR}/external/Dualizer/Dualizer_bin

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/Planarize/Planarize_bin
            ${CPF_DIR}/external/Planarizer/Planarize_bin

            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/extern/gptoolbox-master/mex/fast_sparse.mexa64
            ${CPF_DIR}/external/fast_sparse/fast_sparse.mexa64
            )
  endif()
endif()
