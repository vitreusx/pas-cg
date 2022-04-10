#include <fstream>
#include <ioxx/xyaml/node.h>

namespace ioxx::xyaml {

YAML::Emitter &operator<<(YAML::Emitter &out, node const &node) {
  out << static_cast<YAML::Node const &>(node);
  return out;
}

std::filesystem::path node::abs_path(std::filesystem::path rel_path) const {
  if (loc.has_value())
    return loc.value().parent_path() / rel_path;
  else
    return rel_path;
}

std::filesystem::path node::rel_path(std::filesystem::path abs_path) const {
  if (loc.has_value())
    return std::filesystem::relative(abs_path, loc.value().parent_path());
  else
    return abs_path;
}

node node::import(const std::filesystem::path &path) {
  auto res = node(YAML::LoadFile(path));
  res.loc = path;
  return res;
}

node node::new_file(const std::filesystem::path &path) {
  auto res = node();
  res.loc = path;
  return res;
}

void node::save() const {
  if (loc.has_value()) {
    std::filesystem::create_directories(loc.value().parent_path());
    auto file = std::ofstream(loc.value());
    file << *this;
  }
}

proxy proxy::saving_to(const node &ref) { return proxy(ref, proxy_mode::SAVE); }

proxy proxy::loading_from(const node &ref) {
  return proxy(ref, proxy_mode::LOAD);
}

proxy::proxy(const node &ref, proxy_mode mode) : node(ref), mode(mode) {}

node node::child(const YAML::Node &node) const {
  auto res = ioxx::xyaml::node(node);
  res.loc = loc;
  return res;
}

proxy proxy::child(const YAML::Node &node) const {
  return proxy(this->node::child(node), mode);
}

void node::merge(const node &other) {
  if (other.IsScalar() || !IsMap()) {
    *this = other;
  } else {
    for (auto entry : other) {
      if (entry.first.IsScalar()) {
        auto key = entry.first.Scalar();
        static_cast<node &>(*this)[key].merge(other[key]);
      }
    }
  }
}

node node::clone() const {
  auto clone_ = node(YAML::Clone(*this));
  clone_.loc = loc;
  if (children) {
    clone_.children = std::make_shared<std::unordered_map<std::string, node>>();
    for (auto const &[key, value] : *children) {
      clone_.children->insert(std::make_pair(key, value.clone()));
    }
  }
  return clone_;
}

YAML::Node node::flatten() const {
  auto clone_ = clone();
  if (children) {
    for (auto const &[key, value] : *children) {
      clone_.YAML::Node::operator[](key) = value.flatten();
    }
  }
  return static_cast<YAML::Node const &>(clone_);
}
} // namespace ioxx::xyaml
