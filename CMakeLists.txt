cmake_minimum_required(VERSION 3.5)
project(BAM C)

set(CMAKE_C_STANDARD 11)

set (CMAKE_C_FLAGS "-Wall -pedantic -Werror -o4")
add_executable(BAM src/main.c src/bam.c src/include/bam.h)
target_link_libraries(BAM -lm -lgsl -lgslcblas)
