
/*
    gaussian material class definition (Mandy)
 */


// materials/gaussian.cpp*
#include "materials/gaussian.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "interaction.h"
#include "paramset.h"

namespace pbrt {

// GaussianMaterial Definitions
GaussianMaterial::~GaussianMaterial(){
  delete gm;
}
void GaussianMaterial::ComputeScatteringFunctions(SurfaceInteraction *si,
                                              MemoryArena &arena,
                                              TransportMode mode,
                                              bool allowMultipleLobes) const {
    // Perform bump mapping with _bumpMap_, if present
    if (bumpMap) Bump(bumpMap, si);
    Float e = eta->Evaluate(*si);

    Spectrum op = opacity->Evaluate(*si).Clamp();
    Spectrum t = (-op + Spectrum(1.f)).Clamp();   // default t is 0
    if (!t.IsBlack()) {
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, 1.f);
        BxDF *tr = ARENA_ALLOC(arena, SpecularTransmission)(t, 1.f, 1.f, mode);
        si->bsdf->Add(tr);
    } else
        si->bsdf = ARENA_ALLOC(arena, BSDF)(*si, e);

    Spectrum kd = op * Kd->Evaluate(*si).Clamp();
    if (!kd.IsBlack()) {
        BxDF *diff = ARENA_ALLOC(arena, LambertianReflection)(kd);
        si->bsdf->Add(diff);
    }

    Spectrum ks = op * Ks->Evaluate(*si).Clamp();
    if (!ks.IsBlack()) {
        //Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, e);
        Float roughu, roughv;
        if (roughnessu)
            roughu = roughnessu->Evaluate(*si);
        else
            roughu = roughness->Evaluate(*si);
        if (roughnessv)
            roughv = roughnessv->Evaluate(*si);
        else
            roughv = roughu;

        MicrofacetDistribution *distrib =
          ARENA_ALLOC(arena, BeckmannDistribution)(roughu, roughv, false);
        // use Gaussian BSDF (Mandy)
        BxDF *spec = ARENA_ALLOC(arena, GaussianBSDF)(ks,gm,distrib);
        //std::cout<<"add gaussian bsdf" << "\n";
        //fflush(stdout);
        si->bsdf->Add(spec);
    }

    Spectrum kr = op * Kr->Evaluate(*si).Clamp();
    if (!kr.IsBlack()) {
        Fresnel *fresnel = ARENA_ALLOC(arena, FresnelDielectric)(1.f, e);
        si->bsdf->Add(ARENA_ALLOC(arena, SpecularReflection)(kr, fresnel));
    }

    Spectrum kt = op * Kt->Evaluate(*si).Clamp();
    if (!kt.IsBlack())
        si->bsdf->Add(
            ARENA_ALLOC(arena, SpecularTransmission)(kt, 1.f, e, mode));
}

  GaussianMaterial *CreateGaussianMaterial(const TextureParams &mp) {
    std::shared_ptr<Texture<Spectrum>> Kd =
        mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    std::shared_ptr<Texture<Spectrum>> Ks =
        mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    std::shared_ptr<Texture<Spectrum>> Kr =
        mp.GetSpectrumTexture("Kr", Spectrum(0.f));
    std::shared_ptr<Texture<Spectrum>> Kt =
        mp.GetSpectrumTexture("Kt", Spectrum(0.f));
    std::shared_ptr<Texture<Float>> roughness =
        mp.GetFloatTexture("roughness", .1f);
    std::shared_ptr<Texture<Float>> uroughness =
        mp.GetFloatTextureOrNull("uroughness");
    std::shared_ptr<Texture<Float>> vroughness =
        mp.GetFloatTextureOrNull("vroughness");
    std::shared_ptr<Texture<Float>> eta = mp.GetFloatTextureOrNull("eta");
    if (!eta) eta = mp.GetFloatTexture("index", 1.5f);
    std::shared_ptr<Texture<Spectrum>> opacity =
        mp.GetSpectrumTexture("opacity", 1.f);
    std::shared_ptr<Texture<Float>> bumpMap =
        mp.GetFloatTextureOrNull("bumpmap");
    bool remapRoughness = mp.FindBool("remaproughness", true);

    bool reflectdata = mp.FindBool("reflectdata", true);
    int numgaussian = mp.FindInt("numgaussian", 50);
    float extf = mp.FindFloat("extfactor",1.0);

    // isotropic roughu = roughv
    int dim = 3;
    float alpha = 0.5;
    Gaussianmixture *gm = new Gaussianmixture(dim,numgaussian,alpha,extf);
    return new GaussianMaterial(Kd, Ks, Kr, Kt, roughness, uroughness, vroughness,
                                opacity, eta, bumpMap, remapRoughness,gm);
}

}  // namespace pbrt
