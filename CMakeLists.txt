cmake_minimum_required(VERSION 3.22)

project(haplokit LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)
enable_testing()

include(cmake/htslib.cmake)
haplokit_configure_htslib()

add_library(haplokit_core
    src/cpp/selector.cpp
    src/cpp/vcf_reader.cpp
    src/cpp/view_backend.cpp
    src/cpp/gff_annotator.cpp
    deps/gffsub/src/gff3_parser.cpp
    deps/gffsub/src/gff3_index.cpp
    deps/gffsub/src/gff3_filter.cpp
    deps/gffsub/src/gff3_output.cpp
)

target_include_directories(haplokit_core
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/gffsub/src
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/nlohmann
        ${HAPLOKIT_HTSLIB_INCLUDE_DIRS}
)

target_link_libraries(haplokit_core PUBLIC ${HAPLOKIT_HTSLIB_TARGET} ${HAPLOKIT_HTSLIB_EXTRA_LIBS})
add_dependencies(haplokit_core haplokit_htslib_build)

# Network algorithms library
add_library(haplokit_network STATIC
    src/cpp/network/graph.cpp
    src/cpp/network/hamming_simd.cpp
    src/cpp/network/tcs_algorithm.cpp
    src/cpp/network/msn_algorithm.cpp
    src/cpp/network/mjn_algorithm.cpp
    src/cpp/network/network_backend.cpp
)

target_include_directories(haplokit_network
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/src/cpp/network
)

# Check for AVX2 support
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx2" COMPILER_SUPPORTS_AVX2)
if(COMPILER_SUPPORTS_AVX2)
    target_compile_options(haplokit_network PRIVATE -mavx2)
    target_compile_definitions(haplokit_network PRIVATE HAPLOKIT_USE_SIMD)
endif()

# Check for OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(haplokit_network PUBLIC OpenMP::OpenMP_CXX)
endif()

add_executable(haplokit_cpp
    src/cpp/main.cpp
)

target_link_libraries(haplokit_cpp PRIVATE haplokit_core)

# Network backend executable
add_executable(haplokit_network_backend
    src/cpp/network/main.cpp
)

target_link_libraries(haplokit_network_backend PRIVATE haplokit_network)

if(BUILD_TESTING)
    add_executable(vcf_reader_test
        tests/cpp/vcf_reader_test.cpp
    )

    target_link_libraries(vcf_reader_test PRIVATE haplokit_core)

    add_executable(view_backend_test
        tests/cpp/view_backend_test.cpp
    )

    target_link_libraries(view_backend_test PRIVATE haplokit_core)

    # Network algorithm tests
    add_executable(network_test
        tests/cpp/network_test.cpp
    )

    target_link_libraries(network_test PRIVATE haplokit_network)

    add_test(
        NAME vcf_reader_test
        COMMAND vcf_reader_test
    )

    add_test(
        NAME view_backend_test
        COMMAND view_backend_test
    )

    add_test(
        NAME network_test
        COMMAND network_test
    )

    set_tests_properties(
        vcf_reader_test view_backend_test network_test
        PROPERTIES
            WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    )
endif()
