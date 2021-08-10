cmake_minimum_required(VERSION 3.10)
project(Planarize)

set(CMAKE_CXX_STANDARD 14)

get_directory_property(hasParent PARENT_DIRECTORY)
if(NOT hasParent)
   execute_process(COMMAND git submodule update --init -- ../extern/libigl
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

   execute_process(COMMAND git submodule update --init -- ../extern/directional
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

   execute_process(COMMAND git submodule update --init -- ../extern/ceres-solver
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
                   
   execute_process(COMMAND git submodule update --init -- ../extern/libhedra
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

   execute_process(COMMAND git submodule update --remote
                   WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

   set(BUILD_SHARED_LIBS OFF CACHE BOOL "No shared" FORCE)
   set(BUILD_TESTING OFF CACHE BOOL "No testing targets" FORCE)
   
   if(UNIX AND NOT APPLE)
     set(CGAL_Boost_USE_STATIC_LIBS OFF CACHE BOOL "No static boost libs" FORCE)
   endif()
   
   set(LIBIGL_USE_STATIC_LIBRARY OFF CACHE BOOL "No static libs" FORCE)
   set(LIBIGL_EXPORT_TARGETS ON CACHE BOOL "No igl targets" FORCE)
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
   set(LIBIGL_WITH_MATLAB OFF CACHE BOOL "NO MATLAB" FORCE)
   set(SUITESPARSE OFF CACHE BOOL "NO Eigen SuiteSparse (problems on Win)" FORCE)
   set(CXSPARSE OFF CACHE BOOL "Use CXSPARSE" FORCE)

   set(LIBIGL_WITH_CGAL ON CACHE BOOL "Use CGAL" FORCE)
   set(EIGENSPARSE ON CACHE BOOL "NO Eigen sparse (ceres)" FORCE)
   
   find_package (Eigen3 3.3 REQUIRED NO_MODULE)

  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../extern/libigl ${CMAKE_BINARY_DIR}/extern/libigl)
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../extern/ceres-solver ${CMAKE_BINARY_DIR}/extern/ceres-solver)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../extern/directional/include/ ${CMAKE_BINARY_DIR}/extern/directional)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../extern/libhedra/include/ ${CMAKE_BINARY_DIR}/extern/libhedra)
endif()


# Add your project files
file(GLOB SRCFILES *.cpp)

link_directories(${CMAKE_BINARY_DIR})

add_executable(${PROJECT_NAME}_bin ${SRCFILES})


target_link_libraries(${PROJECT_NAME}_bin igl::core ceres igl::cgal)

