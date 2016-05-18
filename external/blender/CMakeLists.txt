cmake_minimum_required(VERSION 2.6)
project(blender)

include(macros.cmake)

add_subdirectory(intern)

# source after intern and extern to gather all
# internal and external library information first, for test linking
add_subdirectory(source)
