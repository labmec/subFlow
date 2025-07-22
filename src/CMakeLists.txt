#this allow us to define custom FindXXX.cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# This is not the best. Please create (the extense) list of needed files
file(GLOB common_headers *.h json.hpp)
file(GLOB common_sources *.cpp)

# install eigen
if(USING_EIGEN)
    include(cmake/EnableEigen.cmake)
    enable_eigen()
endif() 

# Create the miniMRS library statically!
add_library(miniMRS_library STATIC ${common_headers} ${common_sources})

#enabling TBB library
option(USING_TBB "Whether the TBB library will be linked in" OFF)
if(USING_TBB)
	find_package(Threads)
	find_package(TBB)
	target_compile_definitions(miniMRS_library PRIVATE USING_TBB)
	target_link_libraries(miniMRS_library PUBLIC TBB::tbb)
endif()

# Links the other libraries
if(USING_SUPERLU)    
    target_link_libraries(miniMRS_library PUBLIC ${SUPERLU_LIBRARIES})
endif(USING_SUPERLU)

# Includes the directories
target_include_directories(miniMRS_library PUBLIC ${PZ_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})
if(USING_EIGEN)
    target_include_directories(miniMRS_library PUBLIC ${EIGEN3_INCLUDE_DIR})
endif()
if(USING_SUPERLU)    
    target_include_directories(miniMRS_library PUBLIC ${SUPERLU_INCLUDES})
endif(USING_SUPERLU)

# create folder with input files
set (CONFIG_H "${CMAKE_SOURCE_DIR}/Common/minimrs_config.h")
configure_file ("${CMAKE_SOURCE_DIR}/Common/minimrs_config.h.in" ${CONFIG_H})
