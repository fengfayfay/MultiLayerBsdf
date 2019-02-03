#include "multiscatter.h"
#include "paramset.h"

namespace pbrt {
GaussianMultiScattering* createGaussianMultiScattering(const TextureParams &mp, const std::string& roughness, bool hasTransmission) {
    GaussianMultiScattering *ms = new GaussianMultiScattering;

    bool useMS = mp.FindBool("multiscatter", false);
    bool energyOnly = mp.FindBool("energyonly", false);
    ms->noFresnel = mp.FindBool("noFresnel", true);

    int numgaussian = mp.FindInt("numgaussian", 0);
    float extf = mp.FindFloat("extfactor",1);

    Float alpha = mp.FindFloat(roughness, 0.7f);
    if (useMS) {
        std::cout << "use MultiScatterReflection\n";
        ms->gsReflect = new GaussianScatter(alpha, energyOnly);
        if (hasTransmission) ms->gsTransmit = new GaussianScatter(alpha, energyOnly, true);
    }
    return ms;
}

}//end namespace
