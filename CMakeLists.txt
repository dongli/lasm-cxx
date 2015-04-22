cmake_minimum_required (VERSION 2.8)

project (lasm CXX)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-deprecated-register")

add_subdirectory ("${PROJECT_SOURCE_DIR}/external/geomtk")
include_directories (${GEOMTK_INCLUDE_DIRS})
link_directories (${GEOMTK_LIBRARY_DIRS})

set (source_directories
    "${PROJECT_SOURCE_DIR}/src"
)

foreach (dir ${source_directories})
    include_directories ("${dir}")
    # header files
    file (GLOB tmp1 "${dir}/*.h")
    list (APPEND headers ${tmp1})
    # source files
    aux_source_directory ("${dir}" tmp2)
    list (APPEND sources ${tmp2})
endforeach ()

add_library (lasm ${shared_or_static} ${headers} ${sources})
target_link_libraries (lasm geomtk mlpack)

foreach (dir ${source_directories})
    if (EXISTS "${dir}/test")
        include_directories ("${dir}/test")
        file (GLOB tmp1 "${dir}/test/*.h")
        list (APPEND test_headers ${tmp1})
    endif (EXISTS "${dir}/test")
endforeach ()

include_directories (${gtest_SOURCE_DIR} ${gtest_SOURCE_DIR}/include)
add_executable (test_lasm
    ${test_headers}
    "${PROJECT_SOURCE_DIR}/src/test/test.cpp"
)
target_link_libraries (test_lasm
    gtest gtest_main
    lasm mlpack geomtk
)

add_executable (run_test_cases
    "${PROJECT_SOURCE_DIR}/src/test/run_test_cases.cpp"
)
target_link_libraries (run_test_cases
    lasm
)
