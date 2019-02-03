#ifndef _MULTISCATTER_H
#define _MULTISCATTER_H

#include "gaussianscatter.h"

namespace pbrt{
class GaussianMultiScattering{
public:
    GaussianMultiScattering(): gsReflect(NULL), gsTransmit(NULL) {}
    
    GaussianScatter* gsReflect;
    GaussianScatter* gsTransmit;
    bool noFresnel;
};

GaussianMultiScattering *createGaussianMultiScattering(const TextureParams &mp, const std::string& roughness = "roughness", bool hasTransmission = false);

} //end namespace
#endif
