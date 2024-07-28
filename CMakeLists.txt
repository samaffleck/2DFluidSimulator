cmake_minimum_required(VERSION 3.8)

# Enable Hot Reload for MSVC compilers if supported.
if (POLICY CMP0141)
  cmake_policy(SET CMP0141 NEW)
  set(CMAKE_MSVC_DEBUG_INFORMATION_FORMAT "$<IF:$<AND:$<C_COMPILER_ID:MSVC>,$<CXX_COMPILER_ID:MSVC>>,$<$<CONFIG:Debug,RelWithDebInfo>:EditAndContinue>,$<$<CONFIG:Debug,RelWithDebInfo>:ProgramDatabase>>")
endif()

project(2DFluidSimulator)

# Include directories
#include_directories(include)
include_directories(include libs)

# Add source files
file(GLOB SOURCES "src/*.cpp")

# Add header files
file(GLOB HEADERS "include/2DFluidSimulator/*.h")

# Create static library
add_library(2DFluidSimulator STATIC ${SOURCES} ${HEADERS})

# Set C++ standard after creating the target
if (CMAKE_VERSION VERSION_GREATER 3.12)
  set_property(TARGET 2DFluidSimulator PROPERTY CXX_STANDARD 20)
endif()

# Optionally, if you have tests or other executables, you can add them here
# Example for a test executable:
add_executable(2DFluidSimulatorTest tests/TestMain.cpp)
target_link_libraries(2DFluidSimulatorTest 2DFluidSimulator)
