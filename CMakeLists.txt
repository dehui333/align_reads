cmake_minimum_required(VERSION 3.10)
project (align_reads)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -fPIC")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

execute_process (
    COMMAND which python
    OUTPUT_VARIABLE which_python
)

execute_process (
    COMMAND python3 -c "from sysconfig import get_config_var;print(get_config_var('LIBDIR'), end='/');print(get_config_var('INSTSONAME'), end='')"
    OUTPUT_VARIABLE pylib
)

execute_process (
    COMMAND python3 -c "from sysconfig import get_config_var;print(get_config_var('CONFINCLUDEPY'), end='')"
    OUTPUT_VARIABLE pyinc
)

execute_process (
    COMMAND python3 -c "import numpy as np;print(np.get_include(), end='')"
    OUTPUT_VARIABLE numpyinc
)
message(STATUS "Which python: ${which_python}")
message(STATUS "Python include dir: ${pyinc}")
message(STATUS "Python library path: ${pylib}")
message(STATUS "Numpy include dir: ${numpyinc}")

include_directories(include edlib/include ${pyinc} ${numpyinc})
add_library(align_reads STATIC src/aligner.cpp edlib/src/edlib.cpp)

# bioparser
if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(align_reads bioparser::bioparser ${pylib})

# ram. biosoup and thread_pool contained.
if (NOT TARGET ram)
  add_subdirectory(ram EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(align_reads ram::ram)

# gtest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  tests 
  tests/tests.cpp
)
target_link_libraries(
  tests
  gtest_main
  align_reads
  ${pylib}
)

include(GoogleTest)
gtest_discover_tests(tests)

