#include <iostream>
#include <ioxx/xyaml.h>
using namespace ioxx;

class example {
public:
  std::string a, c;
  xyaml_embed b;

  void connect(xyaml_node_proxy &proxy) {
    proxy["a"] & a;
    proxy["b"] & b;
    b.sub_proxy(proxy)["c"] & c;
  }
};

int main() {
  auto ex = xyaml_node::from_path("data/super.yml").as<example>();

  return 0;
}