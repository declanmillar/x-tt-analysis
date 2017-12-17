cmake_minimum_required (VERSION 3.0 FATAL_ERROR)
project (apsis)

find_package(Boost COMPONENTS system program_options filesystem REQUIRED)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED COMPONENTS RIO Net)
include(${ROOT_USE_FILE})

# find all source files
file(GLOB_RECURSE include include/*.hpp)
file(GLOB_RECURSE src src/*.cpp)
include_directories(include)
include_directories(../delphes/install/include)
find_library(Delphes_LIBRARIES Delphes ../delphes/install/lib)

add_executable(analysis ${src})
target_link_libraries(analysis ${ROOT_LIBRARIES} ${Boost_LIBRARIES} ${Delphes_LIBRARIES})
