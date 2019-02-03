/*
  poly fitted gaussian representation for multi-scattering (Feng)
*/
#if defined(_MSC_VER)
#define NOMINMAX
#pragma once
#endif

#ifndef PBRT_CORE_GAUSSIANSCATTER_H
#define PBRT_CORE_GAUSSIANSCATTER_H

// core/gaussianmixture.h*
#include "pbrt.h"
#include "geometry.h"
#include "stringprint.h"
#include "spectrum.h"
#include <fstream>
#include <vector>

namespace pbrt {

struct Polynomial {

    Polynomial() {}
    Polynomial(Float c0, Float c1, Float c2, Float c3) {
        coef.push_back(c0);
        coef.push_back(c1);
        coef.push_back(c2);
        coef.push_back(c3);
    }

    Float eval(Float x) const {
        Float e = 0;
        for (int i = 0; i< coef.size(); i++) {
            e *= x;
            e += coef[i]; 
        }
        return e;
    }

    Polynomial multiply_x() const {
        Polynomial ipoly;
        if (coef.size() == 0) {
            return ipoly;
        }
        ipoly.coef.resize(coef.size() + 1);
        for (int i = 0; i <coef.size(); i++) {
            ipoly.coef[i] = coef[i]; 
        }
        ipoly.coef[coef.size()] = 0;
        return ipoly;
        
    }
    Polynomial integral() const {
        Polynomial ipoly;
        if (coef.size() == 0) {
            return ipoly;
        }

        ipoly.coef.resize(coef.size() + 1);
        int degree = coef.size() - 1; 
        for (int i = 0; i <= degree; i++) {
            int cdegree = degree - i; 
            ipoly.coef[i] = coef[i]/(cdegree+1); 
        }
        ipoly.coef[coef.size()] = 0;
        return ipoly;
    }

    bool readFromFile(std::ifstream&, int);
    std::vector<Float> coef;
};

  // GaussianScatter class declaration
class GaussianScatter{

  public:
    GaussianScatter();
    GaussianScatter(Float alpha, bool energyOnly = false, bool isTransmission = false);
    ~GaussianScatter();

    bool isEnergyOnly() const { return energyOnly; }

    //hx, hy, wo.z
    Float prob(Float x, Float y, Float z, Float z_i) const {

        if (!validScatter) return 0;
        
        Float energy = penergy.eval(z);
        if (energy < .01 || energyAve < 1e-3) return 0;
        if (energyOnly) {
            Float e_i = penergy.eval(z_i);
            if (e_i < 0) e_i = 0;
            return (energy * e_i) /(M_PI * energyAve);
        }

        Vector2f h(x, y);
        Vector2f mean (pmean[0].eval(z), pmean[1].eval(z));
        Vector2f cov(pcov[0].eval(z), pcov[1].eval(z));
        if (cov[0] < 1e-4 || cov[1] < 1e-4) return 0;
        Vector2f icov(1.0/cov[0], 1.0/cov[1]);

        Vector2f dif = h - mean;
        Vector2f mid(dif[0] * icov[0], dif[1] * icov[1]);
        Float d = Dot(dif, mid);
        Float p = exp(-0.5 * d);

        //gaussian norm factor 2pi^gaussianDimension * determinant of covariance matrix
        Float gnormFactor = 4.0 * M_PI * M_PI * cov[0] * cov[1];
        p *= energy / sqrt(gnormFactor); 
        return p;
    }
    

  private:
    void validateScatter();
    Float alpha;
    bool validScatter;

    Polynomial penergy;
    Polynomial pmean[2];
    Polynomial pcov[2];
    Float energyAve;

    bool energyOnly;
    bool isTransmission;
  };

}// namespace pbrt

#endif  // PBRT_CORE_GAUSSIANSCATTER_H

