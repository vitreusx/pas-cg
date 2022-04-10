#include <ioxx/convert.h>

namespace ioxx {
std::string
convert_impl<std::string, std::string>::operator()(std::string const &s) const {
  return s;
}
} // namespace ioxx
