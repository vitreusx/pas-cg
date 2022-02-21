#pragma once
#include <cg/types/vec3.h>
#include <cg/utils/quantity.h>
#include <ioxx/xyaml.h>

namespace ioxx {
template <> struct xyaml_connection<cg::quantity> {
  void operator()(xyaml_proxy &proxy, cg::quantity &value) const;
};

template <typename Scalar> struct xyaml_connection<cg::vec3<Scalar>> {
  void operator()(xyaml_proxy &proxy, cg::vec3<Scalar> &value) const {
    proxy[0] & value.x();
    proxy[1] & value.y();
    proxy[2] & value.z();
  }
};
} // namespace ioxx