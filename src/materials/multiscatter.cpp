#include "multiscatter.h"
#include "paramset.h"

namespace pbrt {
GaussianMultiScattering* createGaussianMultiScattering(const TextureParams &mp, const std::string& roughness, bool hasTransmission) {
    GaussianMultiScattering *ms = new GaussianMultiScattering;

    bool useMS = mp.FindBool("multiscatter", false);
    bool energyOnly = mp.FindBool("energyonly", false);
    ms->noFresnel = mp.FindBool("noFresnel", true);

    Float alpha = mp.FindFloat(roughness, 0.7f);
    if (useMS) {
        std::cout << "use MultiScatterReflection\n";
        ms->gsReflect = new GaussianScatter(alpha, energyOnly);
        ms->realNVPReflect = new RealNVPScatter();

        if (hasTransmission) ms->gsTransmit = new GaussianScatter(alpha, energyOnly, true);
    }
    return ms;
}

}//end namespace
