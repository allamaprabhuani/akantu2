#include <iostream>
int main() {
  if (true) [[gnu::likely]] {
    std::cout << "has __builtin_expect" << std::endl;
  }
}
