#include "multiscatter.h"
#include "paramset.h"

namespace pbrt {
GaussianMultiScattering* createGaussianMultiScattering(const TextureParams &mp, const std::string& roughness, bool hasTransmission) {
    GaussianMultiScattering *ms = new GaussianMultiScattering;

    bool useMS = mp.FindBool("multiscatter", false);
    bool energyOnly = mp.FindBool("energyonly", false);
    ms->noFresnel = mp.FindBool("noFresnel", false);

    //add find string to nvp model path
    //add find string to nvp model path
    ms->realNVPReflect = NULL;
    ms->gsReflect = NULL;

    Float alpha = mp.FindFloat(roughness, 0.7f);
    bool useNVP = mp.FindBool("usenvp", false);
    int numChannels = mp.FindInt("numChannels", 3);
    
    if (useMS) {
        std::cout << "use MultiScatterReflection\n";
        if (useNVP) {
            std::cout << "use real nvp\n";
            //ms->realNVPReflect = new RealNVPScatter();
            ms->realNVPReflect = new RealNVPScatterSpectrum("exported/gold", numChannels);
        } 
        ms->gsReflect = new GaussianScatter(alpha, energyOnly);
        if (hasTransmission) ms->gsTransmit = new GaussianScatter(alpha, energyOnly, true);
    }
    return ms;
}

}//end namespace
