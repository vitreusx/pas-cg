#pragma once
#include <cg/types/vec3.h>
#include <cg/utils/math.h>

namespace cg {
template <typename U> struct mat3x3 {
  U xx, xy, xz, yx, yy, yz, zx, zy, zz;
  mat3x3(U xx, U xy, U xz, U yx, U yy, U yz, U zx, U zy, U zz)
      : xx{xx}, xy{xy}, xz{xz}, yx{yx}, yy{yy}, yz{yz}, zx{zx}, zy{zy},
        zz{zz} {};
};

template<typename U>
mat3x3<U> operator*(mat3x3<U> const& A, mat3x3<U> const& B) {
  auto xx=A.xx*B.xx+A.xy*B.yx+A.xz*B.zx;
  auto xy=A.xx*B.xy+A.xy*B.yy+A.xz*B.zy;
  auto xz=A.xx*B.xz+A.xy*B.yz+A.xz*B.zz;
  auto yx=A.yx*B.xx+A.yy*B.yx+A.yz*B.zx;
  auto yy=A.yx*B.xy+A.yy*B.yy+A.yz*B.zy;
  auto yz=A.yx*B.xz+A.yy*B.yz+A.yz*B.zz;
  auto zx=A.zx*B.xx+A.zy*B.yx+A.zz*B.zx;
  auto zy=A.zx*B.xy+A.zy*B.yy+A.zz*B.zy;
  auto zz=A.zx*B.xz+A.zy*B.yz+A.zz*B.zz;
  return {xx,xy,xz,yx,yy,yz,zx,zy,zz};
}

template <typename U, typename E>
struct mat_vec_mul_expr : public vec3_expr<mat_vec_mul_expr<U, E>> {
public:
  mat_vec_mul_expr(mat3x3<U> const &A, E const &v) : A{A}, v{v} {};

  decltype(auto) x() const {
    return A.xx * v.x() + A.xy * v.y() + A.xz * v.z();
  }

  decltype(auto) y() const {
    return A.yx * v.x() + A.yy * v.y() + A.yz * v.z();
  }

  decltype(auto) z() const {
    return A.zx * v.x() + A.zy * v.y() + A.zz * v.z();
  }

private:
  mat3x3<U> A;
  E v;
};

template <typename U, typename E>
auto operator*(mat3x3<U> const &A, vec3_expr<E> const &v) {
  return mat_vec_mul_expr<U, E>(A, static_cast<E const &>(v));
}

template <typename E>
mat3x3<double> rot_around(vec3_expr<E> const &a, double angle) {
  auto ct = cos(angle), st = sin(angle);
  auto x = a.x(), y = a.y(), z = a.z();
  auto xx = ct + x * x * (1.0 - ct);
  auto xy = x * y * (1.0 - ct) - z * st;
  auto xz = x * z * (1.0 - ct) + y * st;
  auto yx = y * x * (1.0 - ct) + z * st;
  auto yy = ct + y * y * (1.0 - ct);
  auto yz = y * z * (1.0 - ct) - x * st;
  auto zx = z * x * (1.0 - ct) - y * st;
  auto zy = z * y * (1.0 - ct) + x * st;
  auto zz = ct + z * z * (1.0 - ct);
  return {xx, xy, xz, yx, yy, yz, zx, zy, zz};
}

} // namespace cg