#include "nl/parameters.h"
#include "utils/ioxx_interop.h"
using namespace cg::nl;

void parameters::connect(ioxx::xyaml_proxy &proxy) {
  proxy["pad factor"] & pad_factor;
}