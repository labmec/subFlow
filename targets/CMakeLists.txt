add_executable(test-darcy test-darcy.cpp)
target_link_libraries(test-darcy PRIVATE subFlow_lib)

add_executable(test-transport test-transport.cpp)
target_link_libraries(test-transport PRIVATE subFlow_lib)

add_executable(test-unsaturated-reservoir-vacuum test-unsaturated-reservoir-vacuum.cpp)
target_link_libraries(test-unsaturated-reservoir-vacuum PRIVATE subFlow_lib)

add_executable(test-gravity-segregation test-gravity-segregation.cpp)
target_link_libraries(test-gravity-segregation PRIVATE subFlow_lib)

add_executable(radial-filtration radial-filtration.cpp)
target_link_libraries(radial-filtration PRIVATE subFlow_lib)

add_executable(radial-filtration-3d radial-filtration-3d.cpp)
target_link_libraries(radial-filtration-3d PRIVATE subFlow_lib)

add_executable(three-layers three-layers.cpp)
target_link_libraries(three-layers PRIVATE subFlow_lib)