add_executable(test-darcy test-darcy.cpp)
target_link_libraries(test-darcy PRIVATE subFlow_lib)

add_executable(test-transport test-transport.cpp)
target_link_libraries(test-transport PRIVATE subFlow_lib)

add_executable(test-unsaturated-reservoir-vacuum test-unsaturated-reservoir-vacuum.cpp)
target_link_libraries(test-unsaturated-reservoir-vacuum PRIVATE subFlow_lib)