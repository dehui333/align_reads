cmake_minimum_required(VERSION 3.10)
project(align_reads)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -fPIC")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

execute_process(
  COMMAND which python
  OUTPUT_VARIABLE which_python
)

execute_process(
  COMMAND python3 -c "from sysconfig import get_config_var;print(get_config_var('LIBDIR'), end='')"
  OUTPUT_VARIABLE pylibdir
)

execute_process(
  COMMAND python3 -c "from sysconfig import get_config_var;print(get_config_var('LDVERSION'), end='')"
  OUTPUT_VARIABLE pyversion
)

execute_process(
  COMMAND python3 -c "from sysconfig import get_config_var;print(get_config_var('CONFINCLUDEPY'), end='')"
  OUTPUT_VARIABLE pyinc
)

execute_process(
  COMMAND python3 -c "import numpy as np;print(np.get_include(), end='')"
  OUTPUT_VARIABLE numpyinc
)

message(STATUS "Which python: ${which_python}")

find_package(PythonLibs)

if(EXISTS "${pylibdir}/libpython${pyversion}.so")
  message(STATUS "Using python library: ${pylibdir}/libpython${pyversion}.so")
  message(STATUS "Using python include directory: ${pyinc}")
  set(using_pylib "${pylibdir}/libpython${pyversion}.so")
  set(using_pyinc "${pyinc}")
elseif(EXISTS "${pylibdir}/libpython${pyversion}.a")
  message(STATUS "Using python library: ${pylibdir}/libpython${pyversion}.a")
  message(STATUS "Using python include directory: ${pyinc}")
  set(using_pylib "${pylibdir}/libpython${pyversion}.a")
  set(using_pyinc "${pyinc}")
elseif(EXISTS "${PYTHON_LIBRARIES}" AND NOT IS_DIRECTORY "${PYTHON_LIBRARIES}")
  message(STATUS "Using python library: ${PYTHON_LIBRARIES}")
  message(STATUS "Using python include directory: ${PYTHON_INCLUDE_DIRS}")
  set(using_pylib "${PYTHON_LIBRARIES}")
  set(using_pyinc "${PYTHON_INCLUDE_DIRS}")
else()
  message(FATAL_ERROR "Not able to find python library!")
endif()

# message(STATUS "FOUND ${PYTHON_LIBRARIES}")
# message(STATUS "FOUND2 ${PYTHON_INCLUDE_DIRS}")
# message(STATUS "Python version: ${pyversion}")
# message(STATUS "Python include dir: ${pyinc}")
# message(STATUS "Python library dir: ${pylibdir}")
message(STATUS "Numpy include dir: ${numpyinc}")

# includes
include_directories(include edlib/include ${numpyinc} ${using_pyinc})

# source
add_library(align_reads STATIC edlib/src/edlib.cpp src/Inputs.cpp src/Overlapper.cpp
  src/Generator.cpp
  src/Aligner.cpp
  src/Converter.cpp
  src/AlignmentSegment.cpp
  src/MultiAlignment.cpp
  src/AlignCounter.cpp
  )

# python library
target_link_libraries(align_reads ${using_pylib})

# bioparser
if(NOT TARGET bioparser)
  add_subdirectory(bioparser EXCLUDE_FROM_ALL)
endif()

target_link_libraries(align_reads bioparser::bioparser)

# ram. biosoup and thread_pool contained.
if(NOT TARGET ram)
  add_subdirectory(ram EXCLUDE_FROM_ALL)
endif()

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
  tests/Inputs_tests.cpp
  tests/Overlapper_tests.cpp
  tests/Generator_tests.cpp
  tests/Aligner_test.cpp
  tests/Utilities_test.cpp
  tests/Converter_test.cpp
  tests/AlignmentSegment_test.cpp
  tests/MultiAlignment_test.cpp
  tests/AlignCounter_test.cpp
)
target_link_libraries(
  tests
  gtest_main
  align_reads
)

include(GoogleTest)
gtest_discover_tests(tests)
