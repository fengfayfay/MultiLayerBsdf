/*
  Gaussian mixture for gaussianbsdf (Mandy, Feng)
*/
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_GAUSSIANMIXTURE_H
#define PBRT_CORE_GAUSSIANMIXTURE_H

// core/gaussianmixture.h*
#include "pbrt.h"
#include "geometry.h"
#include "stringprint.h"
#include "spectrum.h"
#include <fstream>
#include <vector>

namespace pbrt {

struct Matrix3x3;

  // Matrix3x3 Declarations
  struct Matrix3x3 {
    // Matrix3x3 Public Methods
    Matrix3x3() {
      m[0][0] = m[1][1] = m[2][2] = 1.f;
      m[0][1] = m[0][2] = m[1][0] = m[1][2] = m[2][0] =
        m[2][1] = 0.f;
    }
    Matrix3x3(Float mat[3][3]);
    Matrix3x3(Float t00, Float t01, Float t02, Float t10, Float t11,
              Float t12, Float t20, Float t21, Float t22);
    Matrix3x3(std::vector<Float> v);
    bool operator==(const Matrix3x3 &m2) const {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if (m[i][j] != m2.m[i][j]) return false;
      return true;
    }
    bool operator!=(const Matrix3x3 &m2) const {
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          if (m[i][j] != m2.m[i][j]) return true;
      return false;
    }
    friend Matrix3x3 Transpose(const Matrix3x3 &);
    void Print(FILE *f) const {
      fprintf(f, "[ ");
      for (int i = 0; i < 3; ++i) {
        fprintf(f, "  [ ");
        for (int j = 0; j < 3; ++j) {
          fprintf(f, "%f", m[i][j]);
          if (j != 2) fprintf(f, ", ");
        }
        fprintf(f, " ]\n");
      }
      fprintf(f, " ] ");
    }
    static Matrix3x3 Mul(const Matrix3x3 &m1, const Matrix3x3 &m2) {
      Matrix3x3 r;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          r.m[i][j] = m1.m[i][0] * m2.m[0][j] + m1.m[i][1] * m2.m[1][j] +
            m1.m[i][2] * m2.m[2][j];
      return r;
    }
    static Vector3f Mul(const Matrix3x3 &m1, const Vector3f&  v) {
      Vector3f r;
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          r[i] = m1.m[i][0] * v[0] + m1.m[i][1] * v[1] + m1.m[i][2] * v[2];
      return r;
    }
    friend Matrix3x3 Inverse(const Matrix3x3 &);
    Float determinant() const{
      Float d = m[0][0]* (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2]*m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1]*m[2][0]);
      return d;
    }

    friend std::ostream &operator<<(std::ostream &os, const Matrix3x3 &m) {
      // clang-format off
      os << StringPrintf("[ [ %f, %f, %f] "
                         "[ %f, %f, %f ] "
                         "[ %f, %f, %f]]",
                         m.m[0][0], m.m[0][1], m.m[0][2],
                         m.m[1][0], m.m[1][1], m.m[1][2],
                         m.m[2][0], m.m[2][1], m.m[2][2]);
      // clang-format on
      return os;
    }


    Float m[3][3];
    Float m_determinant;
  };

  // Gaussianmixture class declaration
  class Gaussianmixture{

  public:
    Gaussianmixture();
    Gaussianmixture(int dim, int num, float alpha, float extf);
    ~Gaussianmixture();
    Float single_gaussian_pdf(Float x, Float y, Float z, int index) const;
    Float prob(Float x, Float y, Float z) const;
    Spectrum brdfcos(const Vector3f &wo, const Vector3f &wi) const;
    void testbrdfcos(int thetanum, int phinum, int munum) const;

  private:
    int dimension;
    int num_gaussian;
    float extfactor;
    std::vector<Float> weights;
    std::vector<Vector3f> means;
    std::vector<Matrix3x3> covars;
    std::vector<Matrix3x3> covars_inverse;

    std::vector<Float> norm_factors;
  };

}// namespace pbrt

#endif  // PBRT_CORE_GAUSSIANMIXTURE_H

