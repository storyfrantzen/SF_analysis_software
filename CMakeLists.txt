cmake_minimum_required(VERSION 3.10)
project(ProcessHipoFiles)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ROOT setup
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

# Environment variable paths
set(CLAS12ROOT $ENV{CLAS12ROOT})
set(HIPO $ENV{HIPO})
set(CCDB_HOME $ENV{CCDB_HOME})
set(RCDB_HOME $ENV{RCDB_HOME})

# Include directories
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${CLAS12ROOT}/include
  ${CLAS12ROOT}/Clas12Banks
  ${CLAS12ROOT}/Clas12Root
  ${HIPO}/include
  ${CCDB_HOME}/include
  ${RCDB_HOME}/cpp/include
)

# Link directories
link_directories(
  ${CLAS12ROOT}/lib
  ${HIPO}/lib
  ${CCDB_HOME}/lib
)

# Collect all .cpp files from src/
file(GLOB SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# List your header files that need dictionaries
set(VARS_HEADERS ${PROJECT_SOURCE_DIR}/include/Vars.h)

# Generate dictionary source file
ROOT_GENERATE_DICTIONARY(G__VarsDict ${VARS_HEADERS} LINKDEF ${PROJECT_SOURCE_DIR}/include/VarsLinkDef.h)

# Add generated dictionary source to your sources
list(APPEND SOURCES ${CMAKE_CURRENT_BINARY_DIR}/G__VarsDict.cxx)

# Executable
add_executable(processHipoFiles ${SOURCES})

# Link libraries
target_link_libraries(processHipoFiles
  ${ROOT_LIBRARIES}
  Clas12Root
  Clas12Banks
  Hipo4
  hipo4
  EG
  pthread
)

