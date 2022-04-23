#pragma once
#include <cg/files/files.h>
#include <cg/types/vec3.h>
#include <cg/utils/quantity.h>

namespace ioxx::xyaml {
template <> struct user_repr<cg::quantity> {
  void load(node const &from, cg::quantity &to) const;
  void save(node &to, cg::quantity const &from) const;
};

template <typename Scalar> struct user_repr<cg::vec3<Scalar>> {
  void link(proxy &proxy, cg::vec3<Scalar> &value) const {
    proxy[0] & value.x();
    proxy[1] & value.y();
    proxy[2] & value.z();
  }
};
} // namespace ioxx::xyaml