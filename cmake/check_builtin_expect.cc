#include <iostream>
int main() {
  if (__builtin_expect(true, 1)) {
    std::cout << "has __builtin_expect" << std::endl;
  }
}
