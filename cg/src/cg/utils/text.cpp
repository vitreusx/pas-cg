#include <algorithm>
#include <cg/utils/text.h>

namespace cg {
std::string lstrip(std::string const &s) {
  auto beg = std::find_if_not(
      s.begin(), s.end(), [](auto c) -> auto { return std::isspace(c); });
  auto end = s.end();
  return s.substr(beg - s.begin(), end - beg);
}

std::string rstrip(std::string const &s) {
  auto beg = s.begin();
  auto end =
      std::find_if_not(
          s.rbegin(), s.rend(), [](auto c) -> auto { return std::isspace(c); })
          .base();
  return s.substr(beg - s.begin(), end - beg);
}

std::string strip(std::string const &s) {
  auto beg = std::find_if_not(
      s.begin(), s.end(), [](auto c) -> auto { return std::isspace(c); });
  auto end =
      std::find_if_not(
          s.rbegin(), s.rend(), [](auto c) -> auto { return std::isspace(c); })
          .base();
  return s.substr(beg - s.begin(), end - beg);
}
} // namespace cg