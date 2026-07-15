cmake_minimum_required(VERSION 3.15...3.26)
project(${SKBUILD_PROJECT_NAME} LANGUAGES CXX)
message(STATUS "Building ${SKBUILD_PROJECT_NAME}")

set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 CONFIG REQUIRED)

find_package(Eigen3)
if(NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found; using github version")
  FetchContent_Declare(
    Eigen3::Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4.0
  )
  FetchContent_MakeAvailable(Eigen3::Eigen)
else()
  message(STATUS "Using system Eigen3")
endif()

add_subdirectory(lib)
pybind11_add_module(_mercury MODULE "src/mercury/binding.cpp" "src/mercury/axial.cpp" "src/mercury/formula1.cpp" "src/mercury/formulav.cpp" "src/mercury/formulaxy.cpp" "src/mercury/lrm.cpp" "src/mercury/rec.cpp")
# uncomment the following line if LRF fitting takes a lot of time
# target_compile_options(_mercury PRIVATE -O2)
target_include_directories(_mercury PRIVATE include)
target_link_libraries(_mercury PRIVATE Mercury)
install(TARGETS _mercury DESTINATION ${SKBUILD_PROJECT_NAME})
