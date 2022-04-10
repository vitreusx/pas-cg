#pragma once
#include "node.h"
#include <optional>

namespace ioxx::xyaml {
template <typename T> struct user_repr<std::optional<T>> {
  void load(node const &from, std::optional<T> &to) const {
    if (from)
      to = from.as<T>();
  }

  void save(node &to, std::optional<T> const &from) const {
    if (from.has_value())
      to << from.value();
  }
};

template <> struct user_repr<std::filesystem::path> {
  void load(node const &from, std::filesystem::path &to) const;
  void save(node &to, std::filesystem::path const &from) const;
};

} // namespace ioxx::xyaml