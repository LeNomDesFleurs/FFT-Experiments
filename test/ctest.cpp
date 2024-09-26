// #include "../src/Tools.hpp"

#include <gtest/gtest.h>

#include <string>
std::string greet(std::string name) { return "Hello " + name + "!"; }
//==============================================================================

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  // Expect two strings not to be equal.
  EXPECT_STRNE("hello", "world");
  // Expect equality.
  EXPECT_EQ(7 * 6, 42);
}

TEST(GreetTest, GreetWithName) {
  EXPECT_EQ(greet("Alice"), "Hello Alice!");
  EXPECT_EQ(greet("Alan Turing"), "Hello Alan Turing!");
}

// int main() { 
//     float one = 1;
//   float two = 2;
//   float three = two + one;
//   float clip = noi::Tools::clipValue(three, 0.f, 2.f);

//   if (clip == 2.f){
//     return 0;
//   } else {
//     return 1;
//   }
// }
