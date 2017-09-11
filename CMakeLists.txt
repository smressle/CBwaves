cmake_minimum_required(VERSION 3.5)

project(CBwaves LANGUAGES CXX)

# Find dependencies
find_package(Threads REQUIRED)
find_package(Boost REQUIRED)

add_executable(${PROJECT_NAME} src/cbwaves.cxx
                               src/ConfigFile.cxx)

target_include_directories(${PROJECT_NAME} PRIVATE include
                                                   ${PROJECT_SOURCE_DIR})

target_link_libraries(${PROJECT_NAME} PRIVATE $<$<BOOL:${WIN32}>:Boost::disable_autolinking> # Boost can auto-link when used with MSVC, but relies on relative paths which CMake does not really like
                                              Threads::Threads
                                              Boost::boost)

target_compile_definitions(${PROJECT_NAME} PRIVATE $<$<BOOL:${WIN32}>:_CRT_SECURE_NO_WARNINGS>  # warning C4996: 'strcpy': This function or variable may be unsafe. Consider using strcpy_s instead.
                                                   $<$<BOOL:${WIN32}>:_SCL_SECURE_NO_WARNINGS>) # warning C4996: 'std::copy::_Unchecked_iterators::_Deprecate': Call to 'std::copy' with parameters that may be unsafe

set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 14
                                                 CXX_STANDARD_REQUIRED ON)