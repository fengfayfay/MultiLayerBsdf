#ifndef _MULTISCATTER_H
#define _MULTISCATTER_H

#include "gaussianscatter.h"
#include "mlbrdf.h"
#include "brdf_lut.h"

namespace pbrt{
class GaussianMultiScattering{
public:
    GaussianMultiScattering(): gsReflect(NULL), gsTransmit(NULL), realNVPReflect(NULL), alpha(0.5) {}
   
    Float alpha; 
    GaussianScatter* gsReflect;
    GaussianScatter* gsTransmit;
    RealNVPScatterSpectrum* realNVPReflect;
    RealNVPScatterSpectrum* realNVPTransmit;
    bool noFresnel;
    bool useBeckmann;
    BrdfLUT *brdfLutAll;
    BrdfLUT *brdfLutMs; 
};

GaussianMultiScattering *createGaussianMultiScattering(const TextureParams &mp, const std::string& roughness = "roughness", bool hasTransmission = false);

} //end namespace
#endif
