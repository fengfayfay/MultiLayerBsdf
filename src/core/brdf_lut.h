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
        fresnel = new Float[tsize * 3];
        for (int i = 0; i < tsize; i++) {
            pdf[i] = 0;
        }
        for (int i = 0; i < tsize * 3; i++) {
            fresnel[i] = 0;
        }
        computeUnits();
    }

    void computeUnits() {
        units[0] = units[1] = 1.0f/(dim-1);
        units[2] = 2.0f * M_PI/(dim-1);
        std::cout <<"units: "<< units << "\n";
    }
    
    void normalize() {
        int pdfvalid = 0;
        for (auto i = 0; i < tsize; i++) {
            auto s = pdf[i];
            if (s > 0) {
                pdfvalid ++;
                for (auto j = 0; j < 3; j++) {
                    if (s > 0) fresnel[i*3 + j] /= s; 
                }
                pdf[i] /= totalSamples;
            }
            
        }
        normed = true;
        std::cout<<"valid pdf: "<< pdfvalid << "\n";
        std::cout<<"total samples: "<< totalSamples << "\n";
    }

    int clamp(Float i) {
        int ind = i;
        if (ind > dim-1) return dim-1;
        if (ind < 0) return 0;
        return i;
    }

    Float computePhi(const Vector3f &wh) {
        auto phi = atan2(wh.y, wh.x);
        if (phi < 0) phi += 2 * M_PI;
        if (phi > 2 * M_PI) phi -= 2 *M_PI;
        if (phi < 0 || phi > 2* M_PI) {
            std::cout<<"phi out of range: " << phi << "\n";
        }
        return phi;
    } 

    uint computeIndex(Float muI, Float muH, Float phiH) {
        if (muI < -1e-3 || muH < -1e-3) {
            std::cout<<"brdf lut expects pos mu\n";
        }
        Float fi = muI/units[0] + 0.5;
        Float fj = muH/units[1] + 0.5;
        Float fk = phiH/units[2] + 0.5;
         
        int i = clamp(fi);
        int j = clamp(fj);
        int k = clamp(fk);
       
        //std::cout << "i j k "<< i << " " << j << " " << k << "\n"; 
        return (i * dim * dim + j * dim + k);
    }   

    uint computeIndex(Float muI, const Vector3f& wh) {
        return computeIndex(muI, wh.z, computePhi(wh));
    }
    uint addSample(Float muI, const Vector3f& wh, Float rgbFresnel[3]) {
        auto i =  computeIndex(muI, wh);
        pdf[i] += 1;
        for (int j = 0; j < 3; j++) {
            fresnel[i*3 + j] += rgbFresnel[j];
        }
        return totalSamples++;
    }

    Float evalPdf(Float muI, const Vector3f& wh) {
        return evalPdf(muI, wh.z, computePhi(wh));
        
    }

    Vector3f evalFresnel(Float muI, const Vector3f& wh) {
        return evalFresnel(muI, wh.z, computePhi(wh));
    }

    Float evalPdf(Float muI, Float muH, Float phiH) {
        return pdf[computeIndex(muI, muH, phiH)];
    }

    Vector3f evalFresnel(Float muI, Float muH, Float phiH) {
        auto i = computeIndex(muI, muH, phiH);
        return Vector3f(fresnel[i*3], fresnel[i*3+1], fresnel[i*3+2]);
    }

    Vector3f evaluateBRDF(Float muI, const Vector3f& wh){
        auto i = computeIndex(muI, wh);
        return Vector3f(fresnel[i * 3], fresnel[i*3 + 1], 
            fresnel[i*3 +2]) * pdf[i];
    }
   
    void save(const std::string &fname);
    void load(const std::string &fname); 

    void setTotalRayCount(uint rayCount) { energyRatio = (Float) totalSamples/rayCount;}
    Float getEnergyRatio() {return energyRatio; }

private:
    Float *pdf;
    Float *fresnel;
    uint dim, rayDepth, totalSamples, tsize;
    Vector3f units;
    Float energyRatio, alpha;
    bool normed;
};

}; //end namespace

#endif
