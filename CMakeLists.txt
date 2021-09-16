cmake_minimum_required(VERSION 3.13)
project("dune-localfunctions" CXX)

# guess dune-common build dir
if(NOT (dune-common_DIR OR dune-common_ROOT OR
        "${CMAKE_PREFIX_PATH}" MATCHES ".*dune-common.*"))
    string(REPLACE  ${CMAKE_PROJECT_NAME} dune-common dune-common_DIR
      ${PROJECT_BINARY_DIR})
endif()

#find dune-common and set the module path
find_package(dune-common REQUIRED)

# start a dune project with information from dune.module
dune_project(dune-localfunctions INTERFACE)

add_subdirectory(doc)
add_subdirectory(dune)

if( DUNE_ENABLE_PYTHONBINDINGS )
  add_subdirectory(python)
  dune_python_install_package(PATH python)
endif()


finalize_dune_project(GENERATE_CONFIG_H_CMAKE)
