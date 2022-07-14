#pragma once
#include "node.h"
#include <cg/types/vec3.h>
#include <cg/utils/quantity.h>
#include <optional>

namespace ioxx::xyaml {
template <typename T> struct user_repr<std::optional<T>> {
  void load(node const &from, std::optional<T> &to) const {
    if (from && (!from.IsScalar() || from.Scalar() != "null"))
      to = from.as<T>();
    else
      to = std::nullopt;
  }

  void save(node &to, std::optional<T> const &from) const {
    if (from.has_value())
      to << from.value();
    else
      to << "null";
  }
};

template <> struct user_repr<std::filesystem::path> {
  void load(node const &from, std::filesystem::path &to) const;
  void save(node &to, std::filesystem::path const &from) const;
};

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