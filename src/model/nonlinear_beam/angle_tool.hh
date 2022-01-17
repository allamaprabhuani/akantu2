#ifndef __ANGLE_TOOL_HPP__
#define __ANGLE_TOOL_HPP__

/**
   This module gives several routines useful to treat angles.
   All is implemented using Eigen, thus trying to reach efficiency by
   vectorisation.
*/

#include "types.hh"
/*
inline Vector18d product_matrix_vector(const Matrix18d &M, const Vector18d &V){
  Vector18d V1;
  for (UInt i=0;i<18;++i)
    for (UInt j=0;i<18;++j)
      V1[i]+=M[i][j] * V[j];
  return V1;
}
*/
inline Matrix3d outer_product(const Vector3d &v1, const Vector3d &v2) {
  /**
      Tensor outer product

      @param v1 vector
      @param v2 vector
      @return Matrix \f$M_{ij} = v^1_i v^2_j\f$
  */

  Matrix3d m;
  for (UInt i = 0; i < 3; ++i)
    for (UInt j = 0; j < 3; ++j)
      m(i, j) = v1[i] * v2[j];
  return m;
}

inline auto skew2vec(const Matrix3d &mat) {
  /**
        From a skew matrix constuct the matching vector

        @param mat A 3x3 skew-symmetric matrix(Matrix)
        @return A 3x1 vector(Matrix)
  */

  auto vec0 = mat(2, 1);
  auto vec1 = mat(0, 2);
  auto vec2 = mat(1, 0);
  Vector3d vec;
  vec << vec0, vec1, vec2;
  return vec;
}

inline Eigen::Matrix3d skew(const Vector3d &vec) {
  /**
      From a vector,
      construct the matching skew matrix

      @param vec Matrix 3x1(vector)
      @return Matrix 3x3
  */
  Eigen::Matrix3d skew;
  skew.row(0) << 0, -vec[2], vec[1];
  skew.row(1) << vec[2], 0, -vec[0];
  skew.row(2) << -vec[1], vec[0], 0;
  return skew;
}

inline Eigen::Matrix3d exp_map(const Vector3d &A) {
  //! Computes the matrix exponential from a vector (representing a rotation)"

  auto theta = A.norm();
  if (theta == 0.) {
    return Eigen::Matrix3d::Identity();
  }
  Eigen::Matrix3d K = skew(A) / theta;

  Eigen::Matrix3d _exp = Eigen::Matrix3d::Identity() + std::sin(theta) * K +
                         (1 - std::cos(theta)) * K * K;
  return _exp;
}

inline auto compute_rotation_angle(const Eigen::Matrix3d &R) {
  //! From a rotation matrix, computes the rotation angle

  auto cos_val = (R.trace() - 1) / 2;
  if (cos_val < -1) {
    cos_val = -1;
  }
  if (cos_val > 1) {
    cos_val = 1;
  }

  auto theta = std::acos(cos_val);
  if (std::isnan(theta)) {
    std::cerr << R << std::endl;
    throw std::runtime_error(
        "in angle calculation theta is not a number: abort");
    // R, (R.trace() - 1) / 2);
  }
  return theta;
}

inline Vector3d log_map(const Eigen::Matrix3d &R) {
  //! Computes the vector matching the rotation passed as argument

  auto theta = compute_rotation_angle(R);
  if (theta < 1e-8) {
    Vector3d res;
    res[0] = -R(1, 2);
    res[1] = R(0, 2);
    res[2] = -R(0, 1);
    return res;
  }
  auto symmetric = double((R - R.transpose()).norm() / R.norm()) < 1e-12;
  if (symmetric) {
    // make it symmetric to avoid numerical errors
    auto _R = 0.5 * (R + R.transpose());
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> s(_R);
    auto eig_vals = s.eigenvalues();
    auto eig_vects = s.eigenvectors();
    for (UInt i = 0; i < eig_vals.size(); ++i) {
      auto val = eig_vals[i];
      auto vect = eig_vects.col(i);
      if (std::abs(val - 1) < 1e-14) {
        return theta * vect / vect.norm();
      }
    }
    throw std::runtime_error(
        "Symmetric cannot log:"); // + std::to_string(R),val);
  }
  // print('R', R)
  auto A = (R - R.transpose()) / 2;
  // print('A', A)
  // print('theta', theta)
  auto _log = theta / std::sin(theta) * A;
  return skew2vec(_log);
}

inline Matrix3d exp_derivative(const Vector3d &B, const Vector3d &Bp) {
  /**
     The purpose of this function is to compute the derivative of a
     rotation function of the kind :

     \f[
      \frac{d}{dt} \Lambda(t)
     \f]

     with \f$\Lambda^T(t) = exp(\theta(t))\f$.

     The implementation is based on Rodrigues's formula.

     @param B the provided rotation vector \f$\theta(t)\f$
     @param Bp the provided rotation vector derivative \f$\frac{d}{dt}
     \theta(t)\f$
  */

  Real f;
  Real h;
  Real fp;
  Real hp;
  Real theta = B.norm();
  if (theta < 1e-12) {
    f = 1;
    h = 0.5;
    fp = 0.;
    hp = 0.;
  } else {
    Real thetap = B.dot(Bp) / theta;
    f = std::sin(theta) / theta;
    h = (1 - std::cos(theta)) / (theta * theta);
    fp = (theta * std::cos(theta) - std::sin(theta)) * thetap / (theta * theta);
    hp = (theta * std::sin(theta) + 2 * std::cos(theta) - 2.) * thetap /
         (theta * theta * theta);
  }
  Matrix3d _B = skew(B);
  Matrix3d _Bp = skew(Bp);
  Matrix3d res = fp * _B + hp * _B * _B + f * _Bp + h * (_Bp * _B + _B * _Bp);
  // print('done exp derive')
  return res;
}

inline Vector3d convected_angle_derivative(const Vector3d &B,
                                           const Vector3d &Bp) {
  /**
     Computes the convected angle derivative .
     This is the closest quantity to an angle variation
     generalized to 3D.

     \f[
        \Lambda^T(t)\frac{d}{dt}\Lambda(t)
     \f]

     with \f$\Lambda^T(t) = exp(\theta(t))\f$.

     The implementation is based on Rodrigues's formula.

     @param B the provided rotation vector \f$\theta(t)\f$
     @param Bp the provided rotation vector derivative
     \f$\frac{d}{dt}\theta(t)\f$

  */

  Matrix3d exp_deriv = exp_map(-B) * exp_derivative(B, Bp);
  Vector3d res = skew2vec(exp_deriv);
  return res;
}

inline auto exp_acceleration(const Vector3d &B, const Vector3d &Bp,
                             const Vector3d &Bpp) {
  auto theta = B.norm();
  auto s = std::sin(theta);
  auto c = std::cos(theta);

  Real f;
  Real h;
  Real fp;
  Real hp;
  Real hpp;
  Real fpp;

  if (theta == 0) {
    f = 1;
    h = 0.5;
    fp = 0;
    hp = 0;
    hpp = 0;
    fpp = 0;
  } else {
    auto thetap = B.dot(Bp) / theta;
    auto thetapp = -std::pow(B.dot(Bp), 2) / std::pow(theta, 3) +
                   (Bp.dot(Bp) + B.dot(Bpp)) / theta;
    f = s / theta;
    h = (1 - c) / std::pow(theta, 2);

    fp = (theta * c - s) * thetap / std::pow(theta, 2);
    hp = (theta * s + 2 * c - 2) * thetap / std::pow(theta, 3);

    fpp = 1 / std::pow(theta, 3) *
          (std::pow(theta, 2) * (thetapp * c - s * std::pow(thetap, 2)) -
           theta * (s * thetapp + 2 * c * std::pow(thetap, 2)) +
           2 * s * std::pow(thetap, 2));

    hpp = 1 / std::pow(theta, 4) *
          (std::pow(theta, 2) * (s * thetapp + c * std::pow(thetap, 2)) +
           2 * theta * ((c - 1) * thetapp - 2 * s * std::pow(thetap, 2)) -
           6 * (c - 1) * std::pow(thetap, 2));
  }
  auto Bhat = skew(B);
  auto Bphat = skew(Bp);
  auto Bpphat = skew(Bpp);
  Eigen::Matrix3d res = fpp * Bhat + hpp * Bhat * Bhat;
  res += 2 * fp * Bphat + 2 * hp * (Bhat * Bphat + Bphat * Bhat);
  res += f * Bpphat + h * (2 * Bphat * Bphat + Bhat * Bpphat + Bpphat * Bhat);
  return res;
}

inline auto angle_acceleration(const Vector3d &B, const Vector3d &Bp,
                               const Vector3d &Bpp) {
  auto L = exp_map(B);
  auto vel = exp_derivative(B, Bp);
  auto acc = exp_acceleration(B, Bp, Bpp);
  auto res = acc * L.transpose() + vel * vel.transpose();
  return skew2vec(res);
}

inline auto compose_rotations(const Vector3d &a1, const Vector3d &a2) {
  /**
      Given two vectors,
    representing rotations,
    this functions
    return the vector representing the composition of the two rotations.*/

  auto a = log_map(exp_map(a1) * exp_map(a2));
  return a;
}

#endif //__ANGLE_TOOL_HPP__
