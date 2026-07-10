cmake_minimum_required(VERSION 3.14)
project(MercuryProject LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(FetchContent)

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

# find_package(ROOT COMPONENTS Minuit2)
if(NOT ROOT_Minuit2_FOUND)
  message(STATUS "Minuit2 not found; using github version")
  FetchContent_Declare(
    Minuit2
    GIT_REPOSITORY https://github.com/root-project/root
    GIT_TAG v6-37-01
    SOURCE_SUBDIR math/minuit2
  )
  FetchContent_MakeAvailable(Minuit2)
  target_compile_definitions(Minuit2Math PUBLIC MATHCORE_STANDALONE)
  set(Minuit_TARGET "Minuit2::Minuit2")

else()
  message(STATUS "Using system ROOT::Minuit2")
  set(Minuit_TARGET "ROOT::Minuit2")
endif()

FetchContent_Declare(
  json11
  GIT_REPOSITORY https://github.com/dropbox/json11.git
  GIT_TAG        master
)
FetchContent_MakeAvailable(json11)

# To to modify the JSON11 library as it is not compatiable with newer gcc compiler
set(json11_cpp "${json11_SOURCE_DIR}/json11.cpp")
file(READ "${json11_cpp}" json11_src)
string(FIND "${json11_src}" "#include <cstdint>" found_include)
if(found_include EQUAL -1)
  string(PREPEND json11_src "#include <cstdint>\n")
  file(WRITE "${json11_cpp}" "${json11_src}")
endif()

file(GLOB_RECURSE MERCURY_CPP_SOURCES CONFIGURE_DEPENDS src/*.cpp)
file(GLOB_RECURSE MERCURY_CPP_EXT_SOURCES CONFIGURE_DEPENDS ext/*.cpp)
set(MERCURY_SOURCES ${MERCURY_CPP_SOURCES} ${MERCURY_CPP_EXT_SOURCES})

add_library(Mercury STATIC ${MERCURY_SOURCES})
target_include_directories(Mercury
  PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include/mercury"
    "${CMAKE_CURRENT_SOURCE_DIR}/ext"
    "${json11_SOURCE_DIR}"
    $<TARGET_PROPERTY:Eigen3::Eigen,INTERFACE_INCLUDE_DIRECTORIES>  # if eigen target available
    "${Minuit2_SOURCE_DIR}/inc"
    "${Minuit2_SOURCE_DIR}/../mathcore/inc"
    "${Minuit2_SOURCE_DIR}/../../core/foundation/inc"
)
target_link_libraries(Mercury PUBLIC json11 Eigen3::Eigen ${Minuit_TARGET})
set_target_properties(Mercury PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(testMercury tests/test.cpp)
target_link_libraries(testMercury PRIVATE Mercury)
target_include_directories(testMercury PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

add_executable(testMercuryXY tests/testxy.cpp)
target_link_libraries(testMercuryXY PRIVATE Mercury)
target_include_directories(testMercuryXY PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")