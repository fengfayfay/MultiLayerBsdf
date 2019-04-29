#ifndef _MULTISCATTER_H
#define _MULTISCATTER_H

#include "gaussianscatter.h"
#include "mlbrdf.h"

namespace pbrt{
class GaussianMultiScattering{
public:
    GaussianMultiScattering(): gsReflect(NULL), gsTransmit(NULL), realNVPReflect(NULL) {}
    
    GaussianScatter* gsReflect;
    GaussianScatter* gsTransmit;
    RealNVPScatterSpectrum* realNVPReflect;
    RealNVPScatterSpectrum* realNVPTransmit;
    bool noFresnel;
};

GaussianMultiScattering *createGaussianMultiScattering(const TextureParams &mp, const std::string& roughness = "roughness", bool hasTransmission = false);

} //end namespace
#endif
