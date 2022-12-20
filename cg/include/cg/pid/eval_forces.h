#pragma once
#include "bundle.h"
#include "lambda.h"
#include <cg/base_forces/lj.h>
#include <cg/base_forces/sink_lj.h>
#include <cg/sbox/pbc.h>
#include <cg/simul/runtime.h>
#include <cg/simul/sched.h>

namespace cg::pid {
class eval_forces : public simul::vect_iter_divisible_mixin<eval_forces> {
  // class eval_forces : public simul::iter_divisible_mixin<eval_forces> {
public:
  lambda bb_plus_lam, bb_minus_lam, ss_lam;
  sink_lj bb_plus_lj, bb_minus_lj;
  vect::vector<sink_lj> ss_ljs;
  real cutoff, active_thr;

public:
  vect::const_view<vec3r> r;
  vect::view<vec3r> F;
  sbox::pbc<real> const *simul_box;
  vect::vector<bundle> const *bundles;
  int const *fast_iter_end;
  real *V, *total_disp;
  vect::const_view<int> prev, next;
  vect::const_view<amino_acid> atype;

public:
  bool is_active(bundle const &bundle) const;

public:
  void iter(int idx) const;
  void fast_iter(int idx) const;
  void vect_iter(int vect_idx) const;
  void fast_vect_iter(int vect_idx) const;
  int size() const;

  class fast_version_t : public simul::divisible {
  public:
    fast_version_t() = default;
    void reset() override;
    void divide(std::vector<simul::task const *> &subtasks,
                int size_hint) override;

  private:
    friend class eval_forces;
    eval_forces const *super = nullptr;
    explicit fast_version_t(eval_forces const *super);

    class slice : public simul::task {
    public:
      explicit slice(eval_forces const *eval, int from, int to);
      void run() const override;

    private:
      eval_forces const *eval;
      int from, to;
    };

    class fast_slice : public simul::task {
    public:
      explicit fast_slice(eval_forces const *eval, int from, int to);
      void run() const override;

    private:
      eval_forces const *eval;
      int from, to;
    };

    class fast_vect_slice : public simul::task {
    public:
      explicit fast_vect_slice(eval_forces const *eval, int from, int to);
      void run() const override;

    private:
      eval_forces const *eval;
      int from, to;
    };

    class vect_slice : public simul::task {
    public:
      explicit vect_slice(eval_forces const *eval, int from, int to);
      void run() const override;

    private:
      eval_forces const *eval;
      int from, to;
    };

    std::vector<slice> slices;
    std::vector<fast_slice> fast_slices;
    std::vector<fast_vect_slice> fast_vect_slices;
    std::vector<vect_slice> vect_slices;
  };

  fast_version_t fast_version() const;

  static inline constexpr int elems_per_vect = vect::VECT_BYTES / sizeof(real);
};
} // namespace cg::pid