#pragma once
#include "../convert.h"
#include <experimental/type_traits>
#include <filesystem>
#include <optional>
#include <unordered_map>
#include <yaml-cpp/yaml.h>

namespace ioxx::xyaml {
class node;
class proxy;
template <typename T> struct user_repr;

class node : public YAML::Node {
public:
  using YAML::Node::Node;

  node(node const &other) = default;
  explicit node(YAML::Node const &other) : YAML::Node(other) {}
  node &operator=(node const &other) = default;

  node(node &&other) = default;
  node &operator=(node &&other) = default;

  static node import(std::filesystem::path const &path);
  static node new_file(std::filesystem::path const &path);

  template <typename Key> node &operator[](Key const &key) {
    if (!children)
      children = std::make_shared<std::unordered_map<std::string, node>>();

    auto key_s = ioxx::convert<std::string>(key);
    if (children->find(key_s) == children->end()) {
      auto yaml_res = this->YAML::Node::operator[](key);
      children->operator[](key_s) = child(yaml_res);
    }
    return children->at(key_s);
  }

  template <typename Key> node operator[](Key const &key) const {
    auto key_s = ioxx::convert<std::string>(key);
    if (children && children->find(key_s) != children->end()) {
      return children->at(key_s);
    } else {
      auto yaml_res = this->YAML::Node::operator[](key);
      return child(yaml_res);
    }
  }

  void merge(node const &other);

  template <typename T> node const &operator>>(T &value) const;

  template <typename T> T as() const {
    T value;
    *this >> value;
    return value;
  }

  template <typename T> node &operator<<(T const &value);

  template <typename T> node &operator=(T const &value) {
    *this << value;
    return *this;
  }

  friend YAML::Emitter &operator<<(YAML::Emitter &out, node const &node);
  void save() const;

  std::filesystem::path abs_path(std::filesystem::path rel_path) const;
  std::filesystem::path rel_path(std::filesystem::path abs_path) const;

  node child(YAML::Node const &node) const;

  node clone() const;
  YAML::Node flatten() const;

public:
  std::optional<std::filesystem::path> loc = std::nullopt;
  mutable std::shared_ptr<std::unordered_map<std::string, node>> children;
};

enum class proxy_mode { LOAD, SAVE };

class proxy : public node {
public:
  proxy_mode mode;
  bool is_loading, is_saving;

  explicit proxy(node const &ref, proxy_mode mode);
  static proxy loading_from(node const &ref);
  static proxy saving_to(node const &ref);

  template <typename Key> proxy operator[](Key const &key) {
    return proxy(this->node::operator[](key), mode);
  }

  template <typename Key> proxy operator[](Key const &key) const {
    return proxy(this->node::operator[](key), mode);
  }

  template <typename T> proxy &operator&(T &value);

  proxy child(YAML::Node const &node) const;
};

template <typename T>
using ext_load_detector = decltype(std::declval<user_repr<T> const &>().load(
    std::declval<node const &>(), std::declval<T &>()));

template <typename T>
inline constexpr bool has_ext_load =
    std::experimental::is_detected_exact_v<void, ext_load_detector, T>;

template <typename T>
using class_load_detector =
    decltype(std::declval<T &>().load(std::declval<node const &>()));

template <typename T>
inline constexpr bool has_class_load =
    std::experimental::is_detected_exact_v<void, class_load_detector, T>;

template <typename T>
using def_load_detector =
    decltype(std::declval<node const &>().YAML::Node::template as<T>());

template <typename T>
inline constexpr bool has_def_load =
    std::experimental::is_detected_v<def_load_detector, T>;

template <typename T>
using ext_save_detector = decltype(std::declval<user_repr<T> const &>().save(
    std::declval<node &>(), std::declval<T const &>()));

template <typename T>
inline constexpr bool has_ext_save =
    std::experimental::is_detected_exact_v<void, ext_save_detector, T>;

template <typename T>
using class_save_detector =
    decltype(std::declval<T const &>().save(std::declval<node &>()));

template <typename T>
inline constexpr bool has_class_save =
    std::experimental::is_detected_exact_v<void, class_save_detector, T>;

template <typename T>
using def_save_detector = decltype(YAML::convert<T>::encode(std::declval<T const&>()));

template <typename T>
inline constexpr bool has_def_save =
    std::experimental::is_detected_v<def_save_detector, T>;

template <typename T>
using ext_link_detector = decltype(std::declval<user_repr<T> const &>().link(
    std::declval<proxy &>(), std::declval<T &>()));

template <typename T>
inline constexpr bool has_ext_link =
    std::experimental::is_detected_exact_v<void, ext_link_detector, T>;

template <typename T>
using class_link_detector =
    decltype(std::declval<T &>().link(std::declval<proxy &>()));

template <typename T>
inline constexpr bool has_class_link =
    std::experimental::is_detected_exact_v<void, class_link_detector, T>;

template <typename T>
inline constexpr bool has_custom_link = has_ext_link<T> || has_class_link<T>;

template <typename T>
inline constexpr bool is_saveable = has_ext_save<T> || has_class_save<T> ||
                                    has_custom_link<T> || has_def_save<T>;

template <typename T>
inline constexpr bool is_loadable = has_ext_load<T> || has_class_load<T> ||
                                    has_custom_link<T> || has_def_load<T>;

template <typename T> struct repr {
  void load(node const &from, T &to) const {
    if constexpr (has_ext_load<T>) {
      user_repr<T>().load(from, to);
    } else if constexpr (has_class_load<T>) {
      to.load(from);
    } else if constexpr (has_custom_link<T>) {
      auto p = proxy::loading_from(from);
      p &to;
    } else if constexpr (has_def_load<T>) {
      to = from.YAML::Node::template as<T>();
    } else {
      throw std::runtime_error("cannot load");
    }
  }

  void save(node &to, T const &from) const {
    if constexpr (has_ext_save<T>) {
      user_repr<T>().save(to, from);
    } else if constexpr (has_class_save<T>) {
      from.save(to);
    } else if constexpr (has_custom_link<T>) {
      auto p = proxy::saving_to(to);
      p &const_cast<T &>(from);
    } else if constexpr (has_def_save<T>) {
      to.YAML::Node::operator=(from);
    } else {
      throw std::runtime_error("cannot save");
    }
  }

  void link(proxy &p, T &x) const {
    if constexpr (has_ext_link<T>) {
      user_repr<T>().link(p, x);
    } else if constexpr (has_class_link<T>) {
      x.link(p);
    } else {
      switch (p.mode) {
      case proxy_mode::LOAD: {
        if constexpr (is_loadable<T>)
          load(p, x);
        else
          throw std::runtime_error("type is not loadable");
        break;
      }
      case proxy_mode::SAVE: {
        if constexpr (is_saveable<T>)
          save(p, x);
        else
          throw std::runtime_error("type is not saveable");
        break;
      }
      }
    }
  }
};

template <typename T> node const &node::operator>>(T &value) const {
  repr<T>().load(*this, value);
  return *this;
}

template <typename T> node &node::operator<<(T const &value) {
  repr<T>().save(*this, value);
  return *this;
}

template <typename T> proxy &proxy::operator&(T &value) {
  repr<T>().link(*this, value);
  return *this;
}

} // namespace ioxx::xyaml