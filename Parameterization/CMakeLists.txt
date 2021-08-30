cmake_minimum_required(VERSION 3.1)
project(Parameterization)

get_directory_property(hasParent PARENT_DIRECTORY)

if(NOT hasParent)
   execute_process(COMMAND git submodule update --init -- ../extern/libigl
                  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

   execute_process(COMMAND git submodule update --init -- ../extern/directional
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
   set(LIBIGL_EXPORT_TARGETS OFF CACHE BOOL "No igl targets" FORCE)
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

   set(LIBIGL_WITH_CGAL ON CACHE BOOL "Use CGAL" FORCE)
   set(LIBIGL_WITH_MATLAB ON CACHE BOOL "Use Matlab" FORCE)

  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../extern/libigl ${CMAKE_BINARY_DIR}/extern/libigl)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../extern/directional/include/ ${CMAKE_BINARY_DIR}/extern/directional)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../extern/libhedra/include/ ${CMAKE_BINARY_DIR}/extern/libhedra)
endif()

if (${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
  #used for matlab to avoid issues
  set (CMAKE_EXE_LINKER_FLAGS -Wl,--copy-dt-needed-entries)
endif()



# Add your project files
file(GLOB SRCFILES *.cpp *.h)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_include_directories(${PROJECT_NAME}_bin PRIVATE ${MATLAB_ROOT}/extern/include)
target_link_libraries(${PROJECT_NAME}_bin PRIVATE igl::core igl::cgal ${CERES_LIBRARIES} igl::matlab ${MATLAB_MX_LIBRARY} ${MATLAB_MEX_LIBRARY} ${MATLAB_MAT_LIBRARY} ${MATLAB_ENG_LIBRARY})
