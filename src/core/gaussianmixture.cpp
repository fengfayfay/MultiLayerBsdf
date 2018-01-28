/*
  gaussianmixture definition (Mandy, Feng)
 */

// core/gaussianmixture.cpp*
#include "gaussianmixture.h"
#include <iostream>
#include <stdarg.h>

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


  Gaussianmixture::Gaussianmixture(int dim, int num, Float alpha, bool reflect): dimension(dim), num_gaussian(num), reflectdata(reflect){
    weights.resize(num_gaussian);
    means.resize(num_gaussian);
    covars.resize(num_gaussian);
    covars_inverse.resize(num_gaussian);
    norm_factors.resize(num_gaussian);

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
        int weightCount = 0;
        while ( getline (weightfile,line) )
          {
            assert(weightCount < num_gaussian);
            weights[weightCount++] = (Float)std::stod(line);
          }
        CHECK_EQ(weightCount, num_gaussian);
        weightfile.close();
      }
    else Error("Unable to open file weights.txt");

    std::cout<<"finished weights\n";

    // means
    std::vector<Vector3f> m;
    std::ifstream meanfile(mf);
    if (meanfile.is_open())
      {
        int lineCount = 0;
        while ( getline (meanfile,line) )
        {
            Float v = (Float)std::stod(line);
            int index = lineCount/num_gaussian;
            means[lineCount%num_gaussian][index] = v;
            lineCount++;
        }
        CHECK_EQ(lineCount, num_gaussian*3);
        meanfile.close();
      }
    else Error("Unable to open file means.txt");
    CHECK_EQ(means.size(), num);
    std::cout<<"finished means\n";

    // covarians
    std::vector<Matrix3x3> c;
    std::vector<Float> c_cur;
    std::ifstream covfile(cf);
    if (covfile.is_open())
      {
        int matCount = 0;
        while ( getline (covfile,line) )
          {
            c_cur.push_back((Float)std::stod(line));
            if (c_cur.size()==dim*dim){
              Matrix3x3 test(c_cur);
              test.m_determinant = test.determinant();
              covars[matCount] = test; 
              covars_inverse[matCount] = Inverse(test);
              matCount++;
              c_cur = {};
            }
          }
        covfile.close();
        CHECK_EQ(matCount, num_gaussian);
      }
    else Error("Unable to open file covars.txt");
    std::cout<<"finished covariant matrices\n";
    fflush(stdout);
   
    //our gaussians is conditioned uniformly over incident angle in (0,.5pi) range 
    Float gaussian_norm_factor = 2.f/M_PI;
    //scale by regular gaussian normalization factor over n dimension 
    for (int i=0; i< dimension; i++) {
        gaussian_norm_factor *= 2.0 * M_PI;
    }

    for (int i = 0; i<num_gaussian; i++) {
        norm_factors[i] = 1.0 / sqrt( gaussian_norm_factor * covars[i].m_determinant);
    }

    // calculate brdf*cos value and write to file for testing
    testbrdfcos(90,400,100);
    std::cout<<"finished output brdfcos\n";
    fflush(stdout);

  }

  Gaussianmixture::~Gaussianmixture(){}

  Float Gaussianmixture::single_gaussian_pdf(Float x, Float y, Float z, int index) const{
    Float p = norm_factors[index];
    std::vector<Float> diff = {x - means[index][0], y-means[index][1], z-means[index][2]};
    std::vector<Float> middle = Matrix3x3::Mul(covars_inverse[index], diff);
    p *= exp(-0.5 * (diff[0]*middle[0] + diff[1]*middle[1] + diff[2]*middle[2]));
    return p;
  }

  Float Gaussianmixture::prob(Float x, Float y, Float z) const{
    Float p = 0.f;
    for (int i = 0; i < num_gaussian; i++){
      p += weights[i] * single_gaussian_pdf(x, y, z, i);
    }
    if (reflectdata) p *= 2;
    return p;
  }

  Spectrum Gaussianmixture:: brdfcos(const Vector3f &wo, const Vector3f &wi) const{
    Float cosThetaO = std::abs(wo.z), cosThetaI = std::abs(wi.z);

    Vector3f wh = wi + wo;

    // handle degenerate cases
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);

    //if (wi.z <=0 || wo.z <= 0) return Spectrum(0.);
    if (wh.z == 0) return Spectrum(0.);

    wh = Normalize(wh);

    // calculate probability in slope domain using fitted GMM
    Float x = wh.x/abs(wh.z);
    Float y = wh.y/abs(wh.z);
    Float zo = acos(abs(wo.z));
    assert(zo >=0 && zo <= M_PI * .5);
    Float p = this->prob(x,y,zo);

    /*
      Float zi = acos(abs(wi.z));
      assert(zi >=0 && zi <= M_PI * .5);
      Float pi = gm->prob(x,y,zi);
      Float recp = fabs(p /cosThetaI - pi/cosThetaO);
      if (recp > 1e-3){
      printf("reciprocity failure\n");
      fflush(stdout);
      }
    */
    
    // probability conditioned on a particular incident angle
    // uniform distributed incident angle prob = 1/(pi/2)
    //p /= 1.f/(M_PI/2);

    Float denom = wi.z + wo.z;
    denom = denom * denom * denom;
    Float J = (wo.z * wi.z + wi.y * wo.y + wi.x * wo.x + 1)/denom;

    return p * J;
  }

  void Gaussianmixture::testbrdfcos(int thetanum, int phinum, int munum) const{
    // wo is incident direction
    Float thetavec[thetanum];
    Float theta_unit = M_PI/2/thetanum;
    for (int i = 0; i < thetanum; ++i){
      thetavec[i] = i * theta_unit;
    }

    // fix wo theta for 60 degree for now
    float angle = 60;
    Vector3f wo = {(float) sin(angle * M_PI/180), 0.f, (float) cos(angle * M_PI/180)};

    // precompute phi mu vec
    Float phivec[phinum];
    Float muvec[munum];
    Float phi_unit = 2 * M_PI / phinum;
    Float mu_unit = 1.0 / munum;

    for (int i = 0; i < phinum; ++i){
      phivec[i] = 2 * M_PI * 0.5 / phinum + phi_unit * i;
    }

    for (int i = 0; i < munum; ++i){
      muvec[i] = 0.5 / munum + mu_unit * i;
    }

    // output file
    std::ofstream output;
    std::ostringstream oss;
    oss << angle << "brdfcos.txt";
    std::string var = oss.str();
    output.open(var);

    for (int i = 0; i < phinum; ++i){
      for (int j = 0; j < munum; ++j){
        // wi is exit direction
        Float phi = phivec[i];
        Float mu = muvec[j];

        Float sintheta = sqrt(1 - mu * mu);
        Vector3f wi = {sintheta * cos(phi), sintheta*sin(phi), mu};
        Spectrum brdfcos = this->brdfcos(wo,wi);

        // write to file
        output << brdfcos[0] <<"\n";
      }
    }
    output.close();
  }




} // namespace pbrt
