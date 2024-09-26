enable_testing()

add_executable(
  hello_test
  ctest.cpp
)
target_link_libraries(
  hello_test
  GTest::gtest_main
)

gtest_discover_tests(hello_test)