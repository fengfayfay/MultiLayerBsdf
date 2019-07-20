#ifndef _BRDFLUT_H
#define _BRDFLUT_H

#include "pbrt.h"
#include "geometry.h"


namespace pbrt{

class BrdfLUT {
public:
    
    BrdfLUT(const std::string& filePrefix, Float alpha, uint rayDepth);
    
    BrdfLUT (Float alpha, uint dim, uint rayDepth) : 
        dim(dim), rayDepth(rayDepth), totalSamples (0), 
        alpha(alpha), normed(false) {
        tsize = dim * dim * dim;
        pdf = new Float[tsize];
        for (int j = 0; j < 3; j++) {
            fresnel[j] = new Float[tsize];
        }
        for (int i = 0; i < tsize; i++) {
            pdf[i] = 0;
            for (int j = 0; j < 3; j++) {
                fresnel[j][i] = 0;
            }
        }
        units[0] = units[1] = 1.0f/(dim-1);
        units[2] = 2.0f * M_PI/(dim-1);
    }
    
    void normalize() {
        for (auto i = 0; i < tsize; i++) {
            auto s = pdf[i];
            if (s > 0) {
                for (auto j = 0; j < 3; j++) {
                    if (s > 0) fresnel[j][i] /= s; 
                }
                s /= totalSamples;
            }
            
        }
        normed = true;
    }

    int clamp(Float ind) {
        if (ind > dim-1) return dim-1;
        if (ind < 0) return 0;
    }

    Float computePhi(const Vector3f &wh) {
        auto phi = atan2(wh.y, wh.x);
        if (phi < 0) phi += 2 * M_PI;
        if (phi > 2 * M_PI) phi -= 2 *M_PI;
    } 

    uint computeIndex(Float muI, Float muH, Float phiH) {
        int i = clamp(muI/units[0]);
        int j = clamp(muH/units[1]);
        int k = clamp(phiH/units[2]);
        return (i * dim * dim + j * dim + k);
    }   

    uint computeIndex(Float muI, const Vector3f& wh) {
        return computeIndex(muI, wh.z, computePhi(wh));
    }
    bool addSample(Float muI, const Vector3f& wh, Float rgbFresnel[3]) {
        auto i =  computeIndex(muI, wh);
        pdf[i] += 1;
        for (int j = 0; j < 3; j++) {
            fresnel[j][i] += rgbFresnel[j];
        }
        totalSamples++;
    }

    Float evalPdf(Float muI, const Vector3f& wh) {
        auto i = computeIndex(muI, wh);
        return pdf[i];
        
    }

    Vector3f evalFresnel(Float muI, const Vector3f& wh) {
        auto i = computeIndex(muI, wh);
        return Vector3f(fresnel[0][i], fresnel[1][i], fresnel[2][i]);
    }

    Float evalPdf(Float muI, Float muH, Float phiH) {
        return pdf[computeIndex(muI, muH, phiH)];
    }

    Vector3f evalFresnel(Float muI, Float muH, Float phiH) {
        auto i = computeIndex(muI, muH, phiH);
        return Vector3f(fresnel[0][i], fresnel[1][i], fresnel[2][i]);
    }

    Vector3f evaluateBRDF(Float muI, const Vector3f& wh){
        auto i = computeIndex(muI, wh);
        return Vector3f(fresnel[0][i], fresnel[1][i], fresnel[2][i]) * pdf[i];
    }
   
    void save(const std::string &fname);
    void load(const std::string &fname); 

    void setTotalRayCount(uint rayCount) { energyRatio = (Float) totalSamples/rayCount;}
    Float getEnergyRatio() {return energyRatio; }

private:
    Float *pdf;
    Float *fresnel[3];
    uint dim, rayDepth, totalSamples, tsize;
    Vector3f units;
    Float energyRatio, alpha;
    bool normed;
};

}; //end namespace

#endif
