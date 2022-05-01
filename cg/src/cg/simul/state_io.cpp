#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/vector.hpp>
#include <cg/simul/state.h>
#include <memory>

namespace boost::serialization {
template <class Archive, typename U>
void serialize(Archive &ar, cg::vec3<U> &v, unsigned int) {
  ar &v.x();
  ar &v.y();
  ar &v.z();
}

template <class Archive, typename U>
void serialize(Archive &ar, cg::box<U> &box, unsigned int) {
  ar &box.cell;
  ar &box.cell_inv;
}

template <class Archive>
void serialize(Archive &ar, cg::sync_data &sync, unsigned int) {
  ar &sync.back();
  ar &sync.side_all();
  ar &sync.side_polar();
  ar &sync.side_hydrophobic();
}

template <typename T> struct data_of {
  data_of() = default;

  explicit data_of(T const &x) {
    auto *data_as_T = reinterpret_cast<T *>(data);
    std::uninitialized_copy_n(&x, 1, data_as_T);
  }

  data_of &operator=(T const &x) {
    auto *data_as_T = reinterpret_cast<T *>(data);
    std::uninitialized_copy_n(&x, 1, data_as_T);
    return *this;
  }

  operator T &() { return *reinterpret_cast<T *>(data); }

  operator T const &() const { return *reinterpret_cast<T const *>(data); }

  alignas(T) unsigned char data[sizeof(T)] = {};
};

template <class Archive, typename T>
void serialize(Archive &ar, data_of<T> &data, unsigned int) {
  ar &data.data;
}

template <typename T> class vector_wrapper {
public:
  explicit vector_wrapper(nitro::vector<T> &v) : v{v} {}

private:
  friend class boost::serialization::access;

  template <class Archive> void load(Archive &ar, unsigned int) {
    std::vector<data_of<T>> stl_v;
    ar >> stl_v;

    v.clear();
    for (int idx = 0; idx < (int)stl_v.size(); ++idx)
      v.emplace_back((T &)stl_v[idx]);
  }

  template <class Archive> void save(Archive &ar, unsigned int) const {
    std::vector<data_of<T>> stl_v;
    stl_v.insert(stl_v.end(), v.begin(), v.end());
    ar &stl_v;
  }

  BOOST_SERIALIZATION_SPLIT_MEMBER()

  nitro::vector<T> &v;
};

#define NV(ar, T, x)                                                           \
  {                                                                            \
    auto wrapper = vector_wrapper<T>(x);                                       \
    ar &wrapper;                                                               \
  }

#define NS(ar, T, x)                                                           \
  {                                                                            \
    auto wrapper = vector_wrapper<nitro::set_node<T>>(x);                      \
    ar &wrapper;                                                               \
  }

template <class Archive>
void serialize(Archive &ar, cg::simul::state &st, unsigned int) {
  ar &st.did_simul_setup;
  ar &st.is_running;

  NV(ar, cg::vec3r, st.r);
  ar &st.box;

  ar &st.t;

  NV(ar, cg::vec3r, st.v);
  NV(ar, cg::vec3sr, st.y0);
  NV(ar, cg::vec3sr, st.y1);
  NV(ar, cg::vec3sr, st.y2);
  NV(ar, cg::vec3sr, st.y3);
  NV(ar, cg::vec3sr, st.y4);
  NV(ar, cg::vec3sr, st.y5);
  ar &st.true_t;

  NS(ar, cg::qa::contact, st.qa_contacts);
  NV(ar, cg::sync_data, st.sync_values);
  NV(ar, bool, st.part_of_ssbond);
  ar &st.num_qa_contacts;

  ar &st.post_equil;
  ar &st.did_post_equil_setup;
}
} // namespace boost::serialization

namespace cg::simul {
std::istream &operator>>(std::istream &is, state &st) {
  using namespace boost::archive;
  binary_iarchive iar(is);
  iar >> st;
  return is;
}

std::ostream &operator<<(std::ostream &os, state const &st) {
  using namespace boost::archive;
  binary_oarchive oar(os);
  oar << st;
  return os;
}
} // namespace cg::simul