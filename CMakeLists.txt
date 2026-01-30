cmake_minimum_required(VERSION 3.14)

project (subFlow LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

# Mac OS X specific code
if (APPLE)
  add_definitions(-DMACOSX)
endif()

# Option to use FastCondensedCompel
option(USE_FAST_CONDENSED_COMPEL "Use FastCondensedCompel" ON)
if(USE_FAST_CONDENSED_COMPEL)
  add_definitions(-DUSE_FAST_CONDENSED_COMPEL)
endif()

# Making pz "seeable" in Xcode
# ----------------
file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})
# ----------------

add_library(subFlow_lib SHARED "")
target_link_libraries(subFlow_lib PUBLIC NeoPZ::pz)
target_include_directories(subFlow_lib PUBLIC ${PZ_INCLUDE_DIRS})
set_target_properties(subFlow_lib PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(src)
add_subdirectory(targets)
