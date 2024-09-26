# to download JUCE
# git submodule update --init --recursive 

cmake_minimum_required(VERSION 3.30)

PROJECT(JUCE_EXAMPLE
    LANGUAGES CXX
    VERSION 0.0.2
)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_OSX_ARCHITECTURES "arm64;x86_64" CACHE STRING "" FORCE)

# for clang-tidy(this enable to find system header files).
if(CMAKE_EXPORT_COMPILE_COMMANDS)
set(CMAKE_CXX_STANDARD_INCLUDE_DIRECTORIES ${CMAKE_CXX_IMPLICIT_INCLUDE_DIRECTORIES})
endif()


# Create a libgtest target to be used as a dependency by test programs

include(GoogleTest)

# Enable JUCE. Do not use find_package to prevent from mix up with one globally installed.
add_subdirectory(lib/JUCE)
add_subdirectory(src)
add_subdirectory(lib/googletest-1.15.2)
add_subdirectory(test)