#include <cg/ckpt/make_checkpoint.h>
#include <cg/utils/text.h>

namespace cg::ckpt {
void make_checkpoint::operator()() const {
  if (*last_t + every <= st->t) {
    std::filesystem::path ckpt_path = format(path_fmt.c_str(), st->t);
    create_directories(ckpt_path.parent_path());
    std::ofstream ckpt_os(ckpt_path);
    ckpt_os << *st;
    *last_t = st->t;
  }
}
} // namespace cg::ckpt