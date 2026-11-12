cmake_minimum_required(VERSION 2.8)

set(CMAKE_BUILD_TYPE Release)
if (UNIX)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -Wextra")
endif(UNIX)
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(redist)

# Boost.MPI already disables the deprecated OpenMPI C++ bindings in its own
# config header.  Keep the equivalent MPICH definition for older installations.
# The Boost.Geometry definition suppresses its informational C++03 notice;
# redist is compiled as C++11 on UNIX.  SYSTEM include directories keep
# third-party header warnings separate from warnings in redist itself.
add_definitions(-DMPICH_SKIP_MPICXX)
add_definitions(-DBOOST_GEOMETRY_DISABLE_DEPRECATED_03_WARNING)

find_package(MPI REQUIRED)

#find_package( CGAL QUIET COMPONENTS  )
#include( ${CGAL_USE_FILE} )

find_package(Boost 1.4 REQUIRED COMPONENTS mpi serialization system timer chrono filesystem)
include_directories(SYSTEM ${MPI_INCLUDE_PATH} ${Boost_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} "redist.cpp")

#TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES} ${CGAL_LIBRARIES} ${MPI_LIBRARIES})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} ${Boost_LIBRARIES} ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(${PROJECT_NAME} PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
