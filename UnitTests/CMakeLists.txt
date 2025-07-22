enable_testing()
#in this file (and ONLY in this file) we will define CATCH_CONFIG_MAIN
add_library(test_library STATIC test_main.cpp)
if (APPLE)
  target_compile_definitions(test_library PUBLIC MACOSX)
endif()

target_link_libraries(test_library PUBLIC NeoPZ::pz iMRS_library Catch2::Catch2)

add_subdirectory(TestDFNIMRS)
add_subdirectory(TestOneFracture)
add_subdirectory(TestOnlyFracIntersec)
add_subdirectory(TestFracIntersect3D)
add_subdirectory(TestNoFrac4Spaces)


# These tests have been discontinued until dfnimrs becomes a class that can be used to call them
# add_subdirectory(TestFlemischCtePressure)

# Mortar spaces have been discontinued in iMRS
# add_subdirectory(TestFracIntersect3DMortar)
