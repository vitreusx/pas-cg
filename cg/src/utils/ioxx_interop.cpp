#include "utils/ioxx_interop.h"
using namespace ioxx;

void xyaml_connection<cg::quantity>::operator()(xyaml_proxy &proxy,
                                                cg::quantity &value) const {
  if (proxy.loading()) {
    std::string repr;
    proxy &repr;
    value = cg::quantity(repr);
  } else {
    auto repr = value.repr();
    proxy &repr;
  }
}
