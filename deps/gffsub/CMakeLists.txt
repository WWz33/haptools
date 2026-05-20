cmake_minimum_required(VERSION 3.10)
project(gffsub VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(gffsub
    src/gffsub.cpp
    src/gff3_index.cpp
    src/gff3_parser.cpp
    src/gff3_filter.cpp
    src/gff3_output.cpp
)

target_include_directories(gffsub PRIVATE
    src
)

target_compile_options(gffsub PRIVATE
    -Wall -O2
)

include(CTest)

if(BUILD_TESTING)
    add_executable(annotation_index_smoke
        tests/annotation_index_smoke.cpp
        src/gff3_index.cpp
        src/gff3_parser.cpp
    )

    target_include_directories(annotation_index_smoke PRIVATE
        src
    )

    target_compile_options(annotation_index_smoke PRIVATE
        -Wall -O2
    )

    add_test(NAME annotation_index_smoke COMMAND annotation_index_smoke)
endif()

# Installation
install(TARGETS gffsub DESTINATION bin)
