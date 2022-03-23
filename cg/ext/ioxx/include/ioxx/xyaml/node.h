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
template <typename T> struct xyaml_conv;

template <typename T>
using load_detector =
    decltype(std::declval<node const &>() >> std::declval<T &>());

template <typename T>
inline constexpr bool is_loadable =
    std::experimental::is_detected_v<load_detector, T>;

template <typename T>
using ext_load_detector = decltype(std::declval<xyaml_conv<T> const &>().load(
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
using save_detector =
    decltype(std::declval<node &>() << std::declval<T const &>());

template <typename T>
inline constexpr bool is_saveable =
    std::experimental::is_detected_v<save_detector, T>;

template <typename T>
using ext_save_detector = decltype(std::declval<xyaml_conv<T> const &>().save(
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
using link_detector = decltype(std::declval<proxy &>() & std::declval<T &>());

template <typename T>
inline constexpr bool is_linkable =
    std::experimental::is_detected_v<link_detector, T>;

template <typename T>
using ext_link_detector = decltype(std::declval<xyaml_conv<T> const &>().link(
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

    auto key_s = convert<std::string>(key);
    if (children->find(key_s) == children->end()) {
      auto yaml_res = this->YAML::Node::operator[](key);
      children->operator[](key_s) = child(yaml_res);
    }
    return children->at(key_s);
  }

  template <typename Key> node operator[](Key const &key) const {
    auto key_s = convert<std::string>(key);
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

template <typename T> node const &node::operator>>(T &value) const {
  using namespace std::experimental;

  if constexpr (has_ext_load<T>) {
    xyaml_conv<T>().load(*this, value);
  } else if constexpr (has_class_load<T>) {
    value.load(*this);
  } else if constexpr (has_custom_link<T>) {
    auto node_proxy = proxy::loading_from(*this);
    node_proxy &value;
  } else {
    value = this->YAML::Node::template as<T>();
  }

  return *this;
}

template <typename T> node &node::operator<<(T const &value) {
  using namespace std::experimental;

  if constexpr (has_ext_save<T>) {
    xyaml_conv<T>().save(*this, value);
  } else if constexpr (has_class_save<T>) {
    value.save(*this);
  } else if constexpr (has_custom_link<T>) {
    auto node_proxy = proxy::saving_to(*this);
    node_proxy &const_cast<T &>(value);
  } else {
    this->YAML::Node::operator=(value);
  }

  return *this;
}

template <typename T> proxy &proxy::operator&(T &value) {
  using namespace std::experimental;
  if constexpr (has_ext_link<T>) {
    xyaml_conv<T>().link(*this, value);
  } else if constexpr (has_class_link<T>) {
    value.link(*this);
  } else {
    switch (mode) {
    case proxy_mode::LOAD:
      *this >> value;
      break;
    case proxy_mode::SAVE:
      *this << value;
      break;
    }
  }

  return *this;
}

} // namespace ioxx::xyaml