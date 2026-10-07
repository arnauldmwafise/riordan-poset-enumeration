cmake_minimum_required(VERSION 3.16)
project(riordan_posets VERSION 1.1.1 LANGUAGES C CXX
        DESCRIPTION "Riordan posets: construction, enumeration and analysis of binary Riordan poset matrices")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(RIORDAN_BUILD_TESTS "Build the test suites" ON)
option(RIORDAN_BUILD_TOOLS "Build the command-line tools (riordan, rcount)" ON)
option(RIORDAN_BUILD_EXAMPLES "Build the examples" ON)

# --- the header-only library ------------------------------------------------------
add_library(riordan INTERFACE)
add_library(riordan::riordan ALIAS riordan)
target_include_directories(riordan INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(riordan INTERFACE cxx_std_20)

if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  set(RIORDAN_WARNINGS -Wall -Wextra -Wpedantic)
endif()

# --- bundled nauty (Apache-2.0), n <= 64, thread-local storage ----------------------
set(NAUTY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/third_party/nauty)
add_library(nauty_w64 STATIC ${NAUTY_DIR}/nauty.c ${NAUTY_DIR}/nautil.c ${NAUTY_DIR}/naugraph.c
            ${NAUTY_DIR}/schreier.c ${NAUTY_DIR}/naurng.c ${NAUTY_DIR}/nautinv.c)
target_compile_definitions(nauty_w64 PUBLIC MAXN=WORDSIZE WORDSIZE=64)
target_include_directories(nauty_w64 PUBLIC ${NAUTY_DIR})

# --- tools --------------------------------------------------------------------------
if(RIORDAN_BUILD_TOOLS)
  find_package(Threads REQUIRED)
  add_executable(riordan_cli tools/riordan_cli.cpp)
  set_target_properties(riordan_cli PROPERTIES OUTPUT_NAME riordan)
  target_link_libraries(riordan_cli PRIVATE riordan)
  target_compile_options(riordan_cli PRIVATE ${RIORDAN_WARNINGS})
  add_executable(rcount tools/rcount.cpp)
  target_link_libraries(rcount PRIVATE riordan nauty_w64 Threads::Threads)
  # Independent verification counter: compiled only against verify/independent, not the library.
  add_executable(independent_rcount verify/independent_rcount.cpp)
  target_include_directories(independent_rcount PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/verify/independent/include)
  target_link_libraries(independent_rcount PRIVATE nauty_w64 Threads::Threads)
  install(TARGETS riordan_cli rcount independent_rcount RUNTIME DESTINATION bin)
endif()

# --- examples -----------------------------------------------------------------------
if(RIORDAN_BUILD_EXAMPLES)
  foreach(ex quickstart sample_listing)
    add_executable(${ex} examples/${ex}.cpp)
    target_link_libraries(${ex} PRIVATE riordan)
    target_compile_options(${ex} PRIVATE ${RIORDAN_WARNINGS})
  endforeach()
endif()

# --- tests (one executable per module, all registered with CTest) -------------------
if(RIORDAN_BUILD_TESTS)
  enable_testing()
  file(GLOB RIORDAN_TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.cpp)
  foreach(src ${RIORDAN_TEST_SOURCES})
    get_filename_component(name ${src} NAME_WE)
    add_executable(${name} ${src})
    target_link_libraries(${name} PRIVATE riordan)
    target_compile_options(${name} PRIVATE ${RIORDAN_WARNINGS})
    add_test(NAME ${name} COMMAND ${name})
  endforeach()
  if(RIORDAN_BUILD_TOOLS)
    add_test(NAME rcount_regression
             COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_rcount.sh $<TARGET_FILE:rcount>)
    add_test(NAME independent_rcount_regression
             COMMAND sh ${CMAKE_CURRENT_SOURCE_DIR}/scripts/check_independent.sh $<TARGET_FILE:independent_rcount>)
  endif()
endif()

install(DIRECTORY include/riordan DESTINATION include)
