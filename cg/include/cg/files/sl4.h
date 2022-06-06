#pragma once
#include "table.h"
#include <algorithm>
#include <memory>
#include <ostream>
#include <unordered_map>
#include <vector>

namespace ioxx::sl4 {
class element {
public:
  virtual ~element() = default;
  virtual void write(std::ostream &os) const = 0;
};

std::ostream &operator<<(std::ostream &os, element const &el);

class div : public element {
public:
  void write(std::ostream &os) const override;

  template <typename T> T &find(int idx) {
    if (idx >= (int)children.size())
      throw std::runtime_error("index " + std::to_string(idx) +
                               " is out of bounds for the children array");
    return *(T *)children[idx].get();
  }

  template <typename T, typename... Args> T &add(Args &&...args) {
    auto ptr = std::make_unique<T>(std::forward<Args>(args)...);
    auto &el_ptr = children.emplace_back(std::move(ptr));
    return *reinterpret_cast<T *>(el_ptr.get());
  }

  template <typename T, typename... Args>
  T &named_add(std::string const &id, Args &&...args) {
    auto &el = add<T>(std::forward<Args>(args)...);
    named_children[id] = &el;
    return el;
  }

  template <typename T> T &find(std::string const &id) {
    if (auto iter = named_children.find(id); iter != named_children.end())
      return *(T *)iter->second;
    else
      throw std::runtime_error("element with name \"" + id + "\" not found");
  }

public:
  std::unordered_map<std::string, element *> named_children;
  std::vector<std::unique_ptr<element>> children;
};

class table : public element {
public:
  table() = default;
  explicit table(ioxx::table::table tab);

  ioxx::table::table *operator->();

  void write(std::ostream &os) const override;

public:
  ioxx::table::table tab;
};

class comment : public element {
public:
  explicit comment(std::string text = "");

  template <typename... Args> explicit comment(Args &&...args) {
    (text += ... += convert<std::string>(args));
    std::for_each(text.begin(), text.end(), [](char &c) {
      c = (char)std::toupper(c);
    });
    text = "#" + text;
  }

  void write(std::ostream &os) const override;

public:
  std::string text;
};

class raw : public element {
public:
  explicit raw(std::string text = "");

  template <typename... Args> explicit raw(Args &&...args) {
    (text += ... += convert<std::string>(args));
  }

  void write(std::ostream &os) const override;

public:
  std::string text;
};
} // namespace ioxx::sl4