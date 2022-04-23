#include <cg/files/xyaml/extras.h>

namespace ioxx::xyaml {

void user_repr<std::filesystem::path>::load(const node &from,
                                            std::filesystem::path &to) const {
  to = from.as<std::string>();
}

void user_repr<std::filesystem::path>::save(
    node &to, const std::filesystem::path &from) const {
  to = from.string();
}

} // namespace ioxx::xyaml