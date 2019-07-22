#include "multiscatter.h"
#include "paramset.h"

namespace pbrt {
MultiScattering* createMultiScattering(const TextureParams &mp, const std::string& roughness, bool hasTransmission) {
    MultiScattering *ms = new MultiScattering;

    bool useMS = mp.FindBool("multiscatter", false);
    bool energyOnly = mp.FindBool("energyonly", false);
    ms->noFresnel = mp.FindBool("noFresnel", false);
    ms ->useBeckmann = mp.FindBool("useBeckman", true);

    //add find string to nvp model path
    //add find string to nvp model path
    ms->realNVPReflect = NULL;
    ms->gsReflect = NULL;
    ms->realNVPTransmit = NULL;
    ms->gsTransmit = NULL;
    ms->brdfLutAll = NULL;
    ms->brdfLutMs = NULL;

    ms->alpha = mp.FindFloat(roughness, 0.7f);
    ms->useLUTMs = mp.FindBool("useLUTMs", false);
    
    std::string brdfPrefix = mp.FindString("brdfPrefix", "None");
    if (brdfPrefix != "None"){
        ms->brdfLutAll = new BrdfLUT(brdfPrefix, ms->alpha, 1);
        ms->brdfLutMs = new BrdfLUT(brdfPrefix, ms->alpha, 2);
    }
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
                ms->gsTransmit = new GaussianScatter(ms->alpha, energyOnly, true); 
            } else {
                ms->gsReflect = new GaussianScatter(ms->alpha, energyOnly);
            }
        }
    }
    return ms;
}

}//end namespace
