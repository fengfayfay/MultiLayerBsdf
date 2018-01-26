/*
  gaussianmixture definition (Mandy)
 */

// core/gaussianmixture.cpp*
#include "gaussianmixture.h"
#include <iostream>

namespace pbrt {
  // Matrix3x3 Method Definitions
  Matrix3x3::Matrix3x3(Float mat[3][3]) { memcpy(m, mat, 9 * sizeof(Float)); }

  Matrix3x3::Matrix3x3(Float t00, Float t01, Float t02, Float t10,
                       Float t11, Float t12, Float t20, Float t21,
                       Float t22) {
    m[0][0] = t00;
    m[0][1] = t01;
    m[0][2] = t02;
    m[1][0] = t10;
    m[1][1] = t11;
    m[1][2] = t12;
    m[2][0] = t20;
    m[2][1] = t21;
    m[2][2] = t22;
  }

  Matrix3x3::Matrix3x3(std::vector<Float> v){
    m[0][0] = v[0];
    m[0][1] = v[1];
    m[0][2] = v[2];
    m[1][0] = v[3];
    m[1][1] = v[4];
    m[1][2] = v[5];
    m[2][0] = v[6];
    m[2][1] = v[7];
    m[2][2] = v[8];
  }

  Matrix3x3 Transpose(const Matrix3x3 &m) {
    return Matrix3x3(m.m[0][0], m.m[1][0], m.m[2][0], m.m[0][1],
                     m.m[1][1], m.m[2][1], m.m[0][2], m.m[1][2],
                     m.m[2][2]);
  }

  Matrix3x3 Inverse(const Matrix3x3 &m) {
    int indxc[3], indxr[3];
    int ipiv[3] = {0, 0, 0};
    Float minv[3][3];
    memcpy(minv, m.m, 3 * 3 * sizeof(Float));
    for (int i = 0; i < 3; i++) {
      int irow = 0, icol = 0;
      Float big = 0.f;
      // Choose pivot
      for (int j = 0; j < 3; j++) {
        if (ipiv[j] != 1) {
          for (int k = 0; k < 3; k++) {
            if (ipiv[k] == 0) {
              if (std::abs(minv[j][k]) >= big) {
                big = Float(std::abs(minv[j][k]));
                irow = j;
                icol = k;
              }
            } else if (ipiv[k] > 1)
              Error("Singular matrix in MatrixInvert");
          }
        }
      }
      ++ipiv[icol];
      // Swap rows _irow_ and _icol_ for pivot
      if (irow != icol) {
        for (int k = 0; k < 3; ++k) std::swap(minv[irow][k], minv[icol][k]);
      }
      indxr[i] = irow;
      indxc[i] = icol;
      if (minv[icol][icol] == 0.f) Error("Singular matrix in MatrixInvert");

      // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
      Float pivinv = 1. / minv[icol][icol];
      minv[icol][icol] = 1.;
      for (int j = 0; j < 3; j++) minv[icol][j] *= pivinv;

      // Subtract this row from others to zero out their columns
      for (int j = 0; j < 3; j++) {
        if (j != icol) {
          Float save = minv[j][icol];
          minv[j][icol] = 0;
          for (int k = 0; k < 3; k++) minv[j][k] -= minv[icol][k] * save;
        }
      }
    }
    // Swap columns to reflect permutation
    for (int j = 2; j >= 0; j--) {
      if (indxr[j] != indxc[j]) {
        for (int k = 0; k < 3; k++)
          std::swap(minv[k][indxr[j]], minv[k][indxc[j]]);
      }
    }
    return Matrix3x3(minv);
  }

  // Gaussianmixture definitions
  Gaussianmixture::Gaussianmixture(){
    dimension = 3;
    num_gaussian = 1;
    weights = {1};
    means = { {0,0,0} };
    covars = {Matrix3x3(1,0,0,0,1,0,0,0,1)};
    reflectdata = false;
  }


  Gaussianmixture::Gaussianmixture(int dim, int num, Float alpha, bool reflect){
    dimension = dim;
    num_gaussian = num;
    reflectdata = reflect;

    std::string wf,mf,cf;
    if (reflect){
      wf = "weights_reflect.txt";
      mf = "means_reflect.txt";
      cf = "covars_reflect.txt";
    }else{
      wf = "weights.txt";
      mf = "means.txt";
      cf = "covars.txt";
    }

    // weights
    std::vector<Float> w;
    std::string line;
    std::ifstream weightfile(wf);
    if (weightfile.is_open())
      {
        while ( getline (weightfile,line) )
          {
            w.push_back((Float)std::stod(line));
          }
        weightfile.close();
      }
    else Error("Unable to open file weights.txt");
    CHECK_EQ(w.size(), num);
    weights = w;

    // means
    std::vector<std::vector<Float>> m;
    std::vector<Float> m_cur;
    std::ifstream meanfile(mf);
    if (meanfile.is_open())
      {
        while ( getline (meanfile,line) )
          {
            m_cur.push_back((Float)std::stod(line));
            if (m_cur.size()==dim){
              m.push_back(m_cur);
              m_cur = {};
            }
          }
        meanfile.close();
      }
    else Error("Unable to open file means.txt");
    CHECK_EQ(m.size(), num);
    means = m;

    // covarians
    std::vector<Matrix3x3> c;
    std::vector<Float> c_cur;
    std::ifstream covfile(cf);
    if (covfile.is_open())
      {
        while ( getline (covfile,line) )
          {
            c_cur.push_back((Float)std::stod(line));
            if (c_cur.size()==dim*dim){
              c.push_back(Matrix3x3(c_cur));
              c_cur = {};
            }
          }
        covfile.close();
      }
    else Error("Unable to open file covars.txt");
    CHECK_EQ(c.size(), num);
    covars = c;

  }

  Gaussianmixture::~Gaussianmixture(){}

  Float Gaussianmixture::single_gaussian_pdf(Float x, Float y, Float z, int index) const
  {
    Float p = 1;
    p *= 1 / sqrt( pow(2*M_PI,dimension) * covars[index].determinant());
    std::vector<Float> diff = {x - means[index][0], y-means[index][1], z-means[index][2]};
    std::vector<Float> middle = Matrix3x3::Mul(Inverse(covars[index]), diff);
    p *= exp(-0.5 * (diff[0]*middle[0] + diff[1]*middle[1] + diff[2]*middle[2]));
    return p;
  }

  Float Gaussianmixture::prob(Float x, Float y, Float z) const
  {
    Float p = 0.f;
    for (int i = 0; i < num_gaussian; i++){
      p += weights[i] * single_gaussian_pdf(x, y, z, i);
    }
    if (reflectdata) p *= 2;
    return p;
  }

} // namespace pbrt
