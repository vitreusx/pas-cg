#pragma once
#include <cstring>
#include <string>

namespace cg {
std::string lstrip(std::string const &s);
std::string rstrip(std::string const &s);
std::string strip(std::string const &s);

template <typename... Args>
std::string format(const char *fmt, Args const &...args) {
  auto req_len = snprintf(nullptr, 0, fmt, args...);
  std::string s(req_len, '\0');
  snprintf(s.data(), req_len + 1, fmt, args...);
  return s;
}

template <typename... Args>
void inplace_format(char *buf, const char *fmt, Args const &...args) {
  auto req_len = snprintf(nullptr, 0, fmt, args...);
  std::string s(req_len, '\0');
  snprintf(s.data(), req_len + 1, fmt, args...);
  memcpy(buf, s.data(), req_len);
}
} // namespace cg