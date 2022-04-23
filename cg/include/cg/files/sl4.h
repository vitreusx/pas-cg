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
  virtual void write(std::ostream &os) = 0;
};

class div : public element {
public:
  void write(std::ostream &os) override;

  element *find(std::string const &id);
  element &operator[](int idx);
  element &operator[](std::string const &id);

  template <typename T, typename... Args> T &add(Args &&...args) {
    auto ptr = std::make_unique<T>(std::forward<Args>(args)...);
    return children.emplace_back(std::move(ptr));
  }

  template <typename T, typename... Args>
  T &add(std::string const &id, Args &&...args) {
    auto &el = add<T>(std::forward<Args>(args)...);
    named_children[id] = &el;
    return el;
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

  void write(std::ostream &os) override;

public:
  ioxx::table::table tab;
};

class comment : public element {
public:
  explicit comment(std::string text = "");

  template <typename... Args> explicit comment(Args &&...args) {
    (text += ... += convert<std::string>(args));
    std::for_each(text.begin(), text.end(),
                  [](char c) -> char { return std::toupper(c); });
  }

  void write(std::ostream &os) override;

public:
  std::string text;
};

class raw : public element {
public:
  explicit raw(std::string text = "");

  template <typename... Args> explicit raw(Args &&...args) {
    (text += ... += convert<std::string>(args));
  }

  void write(std::ostream &os) override;

public:
  std::string text;
};
} // namespace ioxx::sl4