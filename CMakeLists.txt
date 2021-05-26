cmake_minimum_required(VERSION 3.0.0)
project(simple_lk_tsp VERSION 0.1.0)
SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/build)
include_directories(my-help-tools)

include(CTest)
enable_testing()

add_executable(simple_lk_tsp main.cpp TSP.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/Data DESTINATION .)
