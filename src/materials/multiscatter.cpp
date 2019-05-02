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
    ms->realNVPTransmit = NULL;
    ms->gsTransmit = NULL;

    Float alpha = mp.FindFloat(roughness, 0.7f);
    bool useNVP = mp.FindBool("usenvp", false);
    int numChannels = mp.FindInt("numChannels", 3);
    int nfloats = 0;
    //Vector3f energyRatio(0.9575, 0.780155, 0.314);
    Vector3f energyRatio(1.0, 1.0, 1.0);
    energyRatio = mp.FindVector3f("energyRatio", energyRatio);
    
    std::string modelPrefix = mp.FindString("modelPrefix", "exported");
    std::cout<< "path to nvp: "<< modelPrefix << "\n";
    std::string fresnelPrefix = mp.FindString("fresnelPrefix", "None");
    std::cout<< "path to fresnel: "<< fresnelPrefix << "\n";
    
    if (useMS) {
        std::cout << "use MultiScatter\n";
        if (useNVP) {
            std::cout << "use real nvp\n";
            //ms->realNVPReflect = new RealNVPScatter();
            if (hasTransmission) {
                ms->realNVPTransmit = new RealNVPScatterSpectrum(energyRatio, modelPrefix, fresnelPrefix, numChannels);
            } else {
                ms->realNVPReflect = new RealNVPScatterSpectrum(energyRatio, modelPrefix, fresnelPrefix, numChannels);
            }
        } else {
            if (hasTransmission)  {
                ms->gsTransmit = new GaussianScatter(alpha, energyOnly, true); 
            } else {
                ms->gsReflect = new GaussianScatter(alpha, energyOnly);
            }
        }
    }
    return ms;
}

}//end namespace
