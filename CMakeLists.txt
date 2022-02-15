cmake_minimum_required(VERSION 3.10)
project (align_reads)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -fPIC")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(include edlib/include)
add_library(align_reads SHARED src/feature_generator.cpp edlib/src/edlib.cpp)

# bioparser
if (NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif ()
target_link_libraries(align_reads bioparser::bioparser)

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
  tests/basic_tests.cpp
)
target_link_libraries(
  tests
  gtest_main
  align_reads
)

include(GoogleTest)
gtest_discover_tests(tests)