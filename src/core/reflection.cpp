
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

// core/reflection.cpp*
#include "reflection.h"
#include "brdf_lut.h"
#include "multiscatter.h"
#include "spectrum.h"
#include "sampler.h"
#include "sampling.h"
#include "interpolation.h"
#include "scene.h"
#include "interaction.h"
#include "stats.h"
#include <stdarg.h>


namespace pbrt {

// BxDF Utility Functions
Float FrDielectric(Float cosThetaI, Float etaI, Float etaT) {
    cosThetaI = Clamp(cosThetaI, -1, 1);
    // Potentially swap indices of refraction
    bool entering = cosThetaI > 0.f;
    if (!entering) {
        std::swap(etaI, etaT);
        cosThetaI = std::abs(cosThetaI);
    }

    // Compute _cosThetaT_ using Snell's law
    Float sinThetaI = std::sqrt(std::max((Float)0, 1 - cosThetaI * cosThetaI));
    Float sinThetaT = etaI / etaT * sinThetaI;

    // Handle total internal reflection
    if (sinThetaT >= 1) return 1;
    Float cosThetaT = std::sqrt(std::max((Float)0, 1 - sinThetaT * sinThetaT));
    Float Rparl = ((etaT * cosThetaI) - (etaI * cosThetaT)) /
                  ((etaT * cosThetaI) + (etaI * cosThetaT));
    Float Rperp = ((etaI * cosThetaI) - (etaT * cosThetaT)) /
                  ((etaI * cosThetaI) + (etaT * cosThetaT));
    return (Rparl * Rparl + Rperp * Rperp) / 2;
}

// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
Spectrum FrConductor(Float cosThetaI, const Spectrum &etai,
                     const Spectrum &etat, const Spectrum &k) {
    cosThetaI = Clamp(cosThetaI, -1, 1);
    Spectrum eta = etat / etai;
    Spectrum etak = k / etai;

    Float cosThetaI2 = cosThetaI * cosThetaI;
    Float sinThetaI2 = 1. - cosThetaI2;
    Spectrum eta2 = eta * eta;
    Spectrum etak2 = etak * etak;

    Spectrum t0 = eta2 - etak2 - sinThetaI2;
    Spectrum a2plusb2 = Sqrt(t0 * t0 + 4 * eta2 * etak2);
    Spectrum t1 = a2plusb2 + cosThetaI2;
    Spectrum a = Sqrt(0.5f * (a2plusb2 + t0));
    Spectrum t2 = (Float)2 * cosThetaI * a;
    Spectrum Rs = (t1 - t2) / (t1 + t2);

    Spectrum t3 = cosThetaI2 * a2plusb2 + sinThetaI2 * sinThetaI2;
    Spectrum t4 = t2 * sinThetaI2;
    Spectrum Rp = Rs * (t3 - t4) / (t3 + t4);

    return 0.5 * (Rp + Rs);
}

// BxDF Method Definitions
Spectrum ScaledBxDF::f(const Vector3f &wo, const Vector3f &wi) const {
    return scale * bxdf->f(wo, wi);
}

Spectrum ScaledBxDF::Sample_f(const Vector3f &wo, Vector3f *wi,
                              const Point2f &sample, Float *pdf,
                              BxDFType *sampledType) const {
    Spectrum f = bxdf->Sample_f(wo, wi, sample, pdf, sampledType);
    return scale * f;
}

Float ScaledBxDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    return bxdf->Pdf(wo, wi);
}

std::string ScaledBxDF::ToString() const {
    return std::string("[ ScaledBxDF bxdf: ") + bxdf->ToString() +
           std::string(" scale: ") + scale.ToString() + std::string(" ]");
}

Fresnel::~Fresnel() {}
Spectrum FresnelConductor::Evaluate(Float cosThetaI) const {
    return FrConductor(std::abs(cosThetaI), etaI, etaT, k);
}

std::string FresnelConductor::ToString() const {
    return std::string("[ FresnelConductor etaI: ") + etaI.ToString() +
           std::string(" etaT: ") + etaT.ToString() + std::string(" k: ") +
           k.ToString() + std::string(" ]");
}

Spectrum FresnelDielectric::Evaluate(Float cosThetaI) const {
    return FrDielectric(cosThetaI, etaI, etaT);
}

std::string FresnelDielectric::ToString() const {
    return StringPrintf("[ FrenselDielectric etaI: %f etaT: %f ]", etaI, etaT);
}

Spectrum SpecularReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
                                      const Point2f &sample, Float *pdf,
                                      BxDFType *sampledType) const {
    // Compute perfect specular reflection direction
    *wi = Vector3f(-wo.x, -wo.y, wo.z);
    *pdf = 1;
    return fresnel->Evaluate(CosTheta(*wi)) * R / AbsCosTheta(*wi);
}

std::string SpecularReflection::ToString() const {
    return std::string("[ SpecularReflection R: ") + R.ToString() +
           std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
}

Spectrum SpecularTransmission::Sample_f(const Vector3f &wo, Vector3f *wi,
                                        const Point2f &sample, Float *pdf,
                                        BxDFType *sampledType) const {
    // Figure out which $\eta$ is incident and which is transmitted
    bool entering = CosTheta(wo) > 0;
    Float etaI = entering ? etaA : etaB;
    Float etaT = entering ? etaB : etaA;

    // Compute ray direction for specular transmission
    if (!Refract(wo, Faceforward(Normal3f(0, 0, 1), wo), etaI / etaT, wi))
        return 0;
    *pdf = 1;
    Spectrum ft = T * (Spectrum(1.) - fresnel.Evaluate(CosTheta(*wi)));
    // Account for non-symmetry with transmission to different medium
    if (mode == TransportMode::Radiance) ft *= (etaI * etaI) / (etaT * etaT);
    return ft / AbsCosTheta(*wi);
}

std::string SpecularTransmission::ToString() const {
    return std::string("[ SpecularTransmission: T: ") + T.ToString() +
           StringPrintf(" etaA: %f etaB: %f ", etaA, etaB) +
           std::string(" fresnel: ") + fresnel.ToString() +
           std::string(" mode : ") +
           (mode == TransportMode::Radiance ? std::string("RADIANCE")
                                            : std::string("IMPORTANCE")) +
           std::string(" ]");
}

Spectrum LambertianReflection::f(const Vector3f &wo, const Vector3f &wi) const {
    return R * InvPi;
}

std::string LambertianReflection::ToString() const {
    return std::string("[ LambertianReflection R: ") + R.ToString() +
           std::string(" ]");
}

Spectrum LambertianTransmission::f(const Vector3f &wo,
                                   const Vector3f &wi) const {
    return T * InvPi;
}

std::string LambertianTransmission::ToString() const {
    return std::string("[ LambertianTransmission T: ") + T.ToString() +
           std::string(" ]");
}

Spectrum OrenNayar::f(const Vector3f &wo, const Vector3f &wi) const {
    Float sinThetaI = SinTheta(wi);
    Float sinThetaO = SinTheta(wo);
    // Compute cosine term of Oren-Nayar model
    Float maxCos = 0;
    if (sinThetaI > 1e-4 && sinThetaO > 1e-4) {
        Float sinPhiI = SinPhi(wi), cosPhiI = CosPhi(wi);
        Float sinPhiO = SinPhi(wo), cosPhiO = CosPhi(wo);
        Float dCos = cosPhiI * cosPhiO + sinPhiI * sinPhiO;
        maxCos = std::max((Float)0, dCos);
    }

    // Compute sine and tangent terms of Oren-Nayar model
    Float sinAlpha, tanBeta;
    if (AbsCosTheta(wi) > AbsCosTheta(wo)) {
        sinAlpha = sinThetaO;
        tanBeta = sinThetaI / AbsCosTheta(wi);
    } else {
        sinAlpha = sinThetaI;
        tanBeta = sinThetaO / AbsCosTheta(wo);
    }
    return R * InvPi * (A + B * maxCos * sinAlpha * tanBeta);
}

std::string OrenNayar::ToString() const {
    return std::string("[ OrenNayar R: ") + R.ToString() +
           StringPrintf(" A: %f B: %f ]", A, B);
}

Spectrum MicrofacetReflection::f(const Vector3f &wo, const Vector3f &wi) const {
    Float cosThetaO = AbsCosTheta(wo), cosThetaI = AbsCosTheta(wi);
    Vector3f wh = wi + wo;
    // Handle degenerate cases for microfacet reflection
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0.);
    wh = Normalize(wh);
    //flip wh to be on the same side as shading n
    if (wh.z < 0) wh *= -1;

    Spectrum F (1.0);
    if (fresnel) F = fresnel->Evaluate(Dot(wi, wh));
    return R * distribution->D(wh) * distribution->G(wo, wi) * F /
           (4 * cosThetaI * cosThetaO);
}

std::string MicrofacetReflection::ToString() const {
    return std::string("[ MicrofacetReflection R: ") + R.ToString() +
           std::string(" distribution: ") + distribution->ToString() +
           std::string(" fresnel: ") + fresnel->ToString() + std::string(" ]");
}

Spectrum MicrofacetTransmission::f(const Vector3f &wo,
                                   const Vector3f &wi) const {
    if (SameHemisphere(wo, wi)) return 0;  // transmission only

    Float cosThetaO = CosTheta(wo);
    Float cosThetaI = CosTheta(wi);
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0);

    // Compute $\wh$ from $\wo$ and $\wi$ for microfacet transmission
    Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
    Vector3f wh = Normalize(wo + wi * eta);
    if (wh.z < 0) wh = -wh;

    if (Dot(wo, wh) * Dot(wi, wh) > -1e-6) return Spectrum(0);

    Spectrum F = 0;
    if (!noFresnel) F = fresnel.Evaluate(Dot(wo, wh));

    Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
    Float factor = (mode == TransportMode::Radiance) ? (1 / eta) : 1;

    return (Spectrum(1.f) - F) * T *
           std::abs(distribution->D(wh) * distribution->G(wo, wi) * eta * eta *
                    AbsDot(wi, wh) * AbsDot(wo, wh) * factor * factor /
                    (cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
}

std::string MicrofacetTransmission::ToString() const {
    return std::string("[ MicrofacetTransmission T: ") + T.ToString() +
           std::string(" distribution: ") + distribution->ToString() +
           StringPrintf(" etaA: %f etaB: %f", etaA, etaB) +
           std::string(" fresnel: ") + fresnel.ToString() +
           std::string(" mode : ") +
           (mode == TransportMode::Radiance ? std::string("RADIANCE")
                                            : std::string("IMPORTANCE")) +
           std::string(" ]");
}

FresnelBlend::FresnelBlend(const Spectrum &Rd, const Spectrum &Rs,
                           MicrofacetDistribution *distribution)
    : BxDF(BxDFType(BSDF_REFLECTION | BSDF_GLOSSY)),
      Rd(Rd),
      Rs(Rs),
      distribution(distribution) {}
Spectrum FresnelBlend::f(const Vector3f &wo, const Vector3f &wi) const {
    auto pow5 = [](Float v) { return (v * v) * (v * v) * v; };
    Spectrum diffuse = (28.f / (23.f * Pi)) * Rd * (Spectrum(1.f) - Rs) *
                       (1 - pow5(1 - .5f * AbsCosTheta(wi))) *
                       (1 - pow5(1 - .5f * AbsCosTheta(wo)));
    Vector3f wh = wi + wo;
    if (wh.x == 0 && wh.y == 0 && wh.z == 0) return Spectrum(0);
    wh = Normalize(wh);
    Spectrum specular =
        distribution->D(wh) /
        (4 * AbsDot(wi, wh) * std::max(AbsCosTheta(wi), AbsCosTheta(wo))) *
        SchlickFresnel(Dot(wi, wh));
    return diffuse + specular;
}

std::string FresnelBlend::ToString() const {
    return std::string("[ FresnelBlend Rd: ") + Rd.ToString() +
           std::string(" Rs: ") + Rs.ToString() +
           std::string(" distribution: ") + distribution->ToString() +
           std::string(" ]");
}

Spectrum FourierBSDF::f(const Vector3f &wo, const Vector3f &wi) const {
    // Find the zenith angle cosines and azimuth difference angle
    Float muI = CosTheta(-wi), muO = CosTheta(wo);
    Float cosPhi = CosDPhi(-wi, wo);

    // Compute Fourier coefficients $a_k$ for $(\mui, \muo)$

    // Determine offsets and weights for $\mui$ and $\muo$
    int offsetI, offsetO;
    Float weightsI[4], weightsO[4];
    if (!bsdfTable.GetWeightsAndOffset(muI, &offsetI, weightsI) ||
        !bsdfTable.GetWeightsAndOffset(muO, &offsetO, weightsO))
        return Spectrum(0.f);

    // Allocate storage to accumulate _ak_ coefficients
    Float *ak = ALLOCA(Float, bsdfTable.mMax * bsdfTable.nChannels);
    memset(ak, 0, bsdfTable.mMax * bsdfTable.nChannels * sizeof(Float));

    // Accumulate weighted sums of nearby $a_k$ coefficients
    int mMax = 0;
    for (int b = 0; b < 4; ++b) {
        for (int a = 0; a < 4; ++a) {
            // Add contribution of _(a, b)_ to $a_k$ values
            Float weight = weightsI[a] * weightsO[b];
            if (weight != 0) {
                int m;
                const Float *ap = bsdfTable.GetAk(offsetI + a, offsetO + b, &m);
                mMax = std::max(mMax, m);
                for (int c = 0; c < bsdfTable.nChannels; ++c)
                    for (int k = 0; k < m; ++k)
                        ak[c * bsdfTable.mMax + k] += weight * ap[c * m + k];
            }
        }
    }

    // Evaluate Fourier expansion for angle $\phi$
    Float Y = std::max((Float)0, Fourier(ak, mMax, cosPhi));
    Float scale = muI != 0 ? (1 / std::abs(muI)) : (Float)0;

    // Update _scale_ to account for adjoint light transport
    if (mode == TransportMode::Radiance && muI * muO > 0) {
        float eta = muI > 0 ? 1 / bsdfTable.eta : bsdfTable.eta;
        scale *= eta * eta;
    }
    if (bsdfTable.nChannels == 1)
        return Spectrum(Y * scale);
    else {
        // Compute and return RGB colors for tabulated BSDF
        Float R = Fourier(ak + 1 * bsdfTable.mMax, mMax, cosPhi);
        Float B = Fourier(ak + 2 * bsdfTable.mMax, mMax, cosPhi);
        Float G = 1.39829f * Y - 0.100913f * B - 0.297375f * R;
        Float rgb[3] = {R * scale, G * scale, B * scale};
        return Spectrum::FromRGB(rgb).Clamp();
    }
}

std::string FourierBSDF::ToString() const {
    return StringPrintf("[ FourierBSDF eta: %f mMax: %d nChannels: %d nMu: %d ",
                        bsdfTable.eta, bsdfTable.mMax, bsdfTable.nChannels,
                        bsdfTable.nMu) +
           std::string(" mode : ") +
           (mode == TransportMode::Radiance ? std::string("RADIANCE")
                                            : std::string("IMPORTANCE")) +
           std::string(" ]");
}

bool FourierBSDFTable::GetWeightsAndOffset(Float cosTheta, int *offset,
                                           Float weights[4]) const {
    return CatmullRomWeights(nMu, mu, cosTheta, offset, weights);
}

Float computeJacobian(const Vector3f &wo, const Vector3f &wi, Float etaO, Float etaI) {
    Vector3f h = etaI * wi + etaO * wo;

    //if (fabs(h.z) < 1e-6 or fabs(wo.z) < 1e-6) return  0;

    /*
    float denom = wi.z * h.z * h.z;
    denom *= denom;
    if (denom < 1e-6) return 1;
    
    Float dxdxi = etaI * (h.z * wi.z + h.x * wi.x);
    Float dydyi = etaI * (h.z * wi.z + h.y * wi.y);
    Float dxdyi = etaI * (h.x * wi.y);
    Float dydxi = etaI * (h.y * wi.x);
    Float J = (dxdxi * dydyi - dydxi * dxdyi)/denom;
    J *= wi.z;
   
    */

    //denom = (wo(3) + wi(3));
    //jacobian = abs(wo(3)* wi(3)  + wo(2)*wi(2) + wo(1)*wi(1) + 1)/denom; 

    Float denom = h.z;
    denom = denom * denom * denom;
    //if (denom < 1e-6) return 0;
    Float jacobian = fabs((Dot(wo, wi) + 1.0) / denom);
    return jacobian;
}

Spectrum computeMultiScattering(Vector3f &wo, Vector3f &wi, 
    Float etaO, Float etaI, const MultiScattering *ms,
    const MicrofacetDistribution *distribution) {

    if (ms->flip && !ms->useLUTMs  && !ms->useLUTAll) {
        Vector3f tmp = wo;
        wo = wi;
        wi = tmp;
    }
    Float cosThetaI = AbsCosTheta(wi);
    Float cosThetaO = AbsCosTheta(wo);

    //std::cout << "\nwo: "<< wo << " wi: "<< wi <<"\n";

    // handle degenerate cases
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0);

    Float phi = atan2(wo.y, wo.x);
    if (fabs(phi) > 1e-5) {
        phi = Degrees(phi);
        Transform rotation = RotateZ(-phi);
        wo = rotation(wo);
        if (fabs(wo.y) > 1e-3 || fabs(wo.z) < 1e-6) {
            //std::cout << "phi: " << phi << "mu: " << cosThetaO << "\n";
            //std::cout<< "wo: " << wo << "\n";
            //std::cout<< "rotation:" << rotation << "\n";
            //fflush(stdout);
        }
        assert(fabs(wo.y) < 1e-5);
        wi = rotation(wi);
    }


    if (wo.z < 0) {
        wo = -wo;
        wi = -wi;
    }
    Float thetaO = acos(fabs(wo.z));
    if (thetaO < 1e-4) thetaO = 1e-4;
    if (thetaO > .5 * M_PI -1e-5) thetaO =  .5 * M_PI -1e-5;

    
    Vector3f wh = Normalize(etaI * wi + etaO*wo);
    if (fabs(wh.z) < 1e-5) return 0;
    if (fabs(cosThetaO - 1) < 1e-5 and cosThetaI < 1e-5) return 0;
    Float J = computeJacobian(wo, wi, etaO, etaI); 

    Float lutThetaI = cosThetaO;
    if (ms->useLUTAll) {
        //std::cout<<"use lut all\n";
        Float pdf = ms->brdfLutAll->evalPdf(lutThetaI, wh) * J / cosThetaI;
        pdf *= 200 * 200/(M_PI * 2);
        Vector3f brdf = ms->brdfLutAll->evalFresnel(lutThetaI, wh) * pdf;
        return SpectrumFromVector3f(brdf);
    }

    if (ms->useLUTMs) {
        //std::cout<<"use lut ms\n";
        Float pdf = ms->brdfLutMs->evalPdf(lutThetaI, wh) * J/cosThetaI;
        Vector3f brdf = ms->brdfLutAll->evalFresnel(lutThetaI, wh) * pdf;
        return SpectrumFromVector3f(brdf);
    }

    
    // calculate probability in slope domain using fitted GMM
    Float x = wh.x/wh.z;
    Float y = wh.y/wh.z;

    Float sbrdf =  (distribution->D(wh) * distribution->G(wo, wi) /
           (4 * cosThetaI * cosThetaO)) ;

    if (ms->realNVPReflect) {
        //std::cout<<"use nvp\n";
        Float thetaH = acos(wh.z);
        Float phiH = atan2(wh.y, wh.x);
        if (phiH < -M_PI) {
            phiH += 2.0 * M_PI;
        } else {
            if (phiH > M_PI) phiH -= 2.0 * M_PI;
        }
        Spectrum p = ms->realNVPReflect->eval(ms->alpha, thetaO, thetaH, phiH, Vector2f(x, y)); 
        Spectrum f(1.0);
        if (ms->brdfLutAll) {
            f = SpectrumFromVector3f(ms->brdfLutAll->evalFresnel(lutThetaI, wh));
        }
        Spectrum multiS =  p * f *  J  / cosThetaI; 
        return multiS;

    } else {

        //std::cout<<"use gs\n";
        Float zo = fabs(wo.z);
        zo = zo > 1? 1: zo;
        Float zi = fabs(wi.z);
        zi = zi > 1? 1: zi;
        Float p = 0;
        if (ms->gsReflect != NULL) {
            p = ms->gsReflect->prob(x,y, zo, zi);
        }
        Float multi = 0;
        if (ms->gsReflect->isEnergyOnly()) {
            multi = p /cosThetaI;
        } else { 
            multi =  p * J / cosThetaI;
        }

        if (ms->brdfLutMs) {
            auto f = ms->brdfLutMs->evalFresnel(lutThetaI, wh);
            //std::cout << "lut fresnel: " << f << "\n";
            Spectrum fn = SpectrumFromVector3f(f);
            return fn * multi;
        }
        return  multi;
    }
    return 0;
}

////////////////////////////////////////////////////////////////////////////////////////////
//////MultiScatterReflection:  Feng
////////////////////////////////////////////////////////////////////////////////////////////
Spectrum MultiScatterReflection::f(const Vector3f &woO, const Vector3f &wiO) const {
    Spectrum singleScatter (0);
    if (!ms->realNVPReflect || ms->useLUTMs) singleScatter = MicrofacetReflection::f(woO, wiO);

    Vector3f wo = Normalize(woO);
    Vector3f wi = Normalize(wiO);
    Vector3f wh = Normalize(wo + wi);
    
    if (wo.z * wi.z < 0) return 0;

    Spectrum MF(1.0);
    /*
    Spectrum F(1.0);
    if (getFresnel() && !ms->realNVPReflect) {
        if (!brdfLutMs) {
        //F = getFresnel()->Evaluate(Dot(wi, wh));
        F = getFresnel()->Evaluate(0.1);
        //auto F2 = F * F * F;
        auto F2 = F * F * F * F;
        auto F3 = F2 * F2;
        auto F4 = F3 * F3;
        //auto c0 = alpha;
        //auto a2 = 1.f/alpha * alpha;
        //a2 = exp(-a2) * a2;
        //auto a2 = alpha * alpha;
        auto a2 = alpha;
        auto nw = 1.f/(1.0 + a2 + a2 * a2);
        auto c1 = nw;
        auto c2 = a2 * nw;
        auto c3 = a2 * a2 * nw;
         
        //MF =  F2 * 0.8 + F3 * 0.18  +  F4 * .02;
        MF =  F2 * c1 + F3 * c2  +  F4 * c3;
        }
    }
    */

    if (wo.z * wi.z < 0) { 
        return singleScatter; 
    }

    Spectrum multi = computeMultiScattering(wo, wi, 1, 1, ms, distribution);

    if (ms->realNVPReflect || ms->useLUTAll) return multi;
    else return  singleScatter +  MF * multi;
}

////////////////////////////////////////////////////////////////////////////////////////////
//////MultiScatterTransmission:  Feng
////////////////////////////////////////////////////////////////////////////////////////////
Spectrum MultiScatterTransmission::f(const Vector3f &woO, const Vector3f &wiO) const {
    Vector3f wo = Normalize(woO);
    Vector3f wi = Normalize(wiO);
    
    if (SameHemisphere(wo, wi)) return 0;  // transmission only

    Float cosThetaO = CosTheta(wo);
    Float cosThetaI = CosTheta(wi);
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0);

    // Compute $\wh$ from $\wo$ and $\wi$ for microfacet transmission
    bool entering = CosTheta(wo) > 0; 
    Float etaI = entering ? etaA : etaB;
    Float etaT = entering ? etaB : etaA;

    Float eta = entering ? (etaB / etaA) : (etaA / etaB);
    Vector3f wh = Normalize(etaI * wo +  etaT * wi);
    if (wh.z < 0) wh = -wh;
    
    if (Dot(wo, wh) * Dot(wi, wh) > -1e-6) return Spectrum(0);

    //Float FrDielectric(Float cosThetaI, Float etaI, Float etaT) 
    Float F = FrDielectric(Dot(wo, wh), etaA, etaB);
    if (!noFresnel) F =  1.0 - F;
    if (F < 0) F = 0;
   
    /* 
    Float sqrtDenom = Dot(wo, wh) * etaI  + etaT * Dot(wi, wh);
    Spectrum singleScatter =  F  *  T *
           std::abs(distribution->D(wh) * distribution->G(wo, wi) * etaT * etaT*
                    AbsDot(wi, wh) * AbsDot(wo, wh) /
                    (cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));
    
    Spectrum multi(singleScatter);
    */

    //factor term cancel out eta * eta in the single scatter evaluation function
    //original equation from walter paper is FDG * absdot(wi, wh) * absdot(wo, wh) / {dot(wi,n) * dot(wo, n) * sqrtDenom^2}

    Spectrum msBrdf =  F * T *  
        computeMultiScattering(entering? wo:wi, entering? wi: wo, etaI, etaT, ms, 
            distribution);
    Spectrum multi(msBrdf);
    return multi;
}


Spectrum BxDF::Sample_f(const Vector3f &wo, Vector3f *wi, const Point2f &u,
                        Float *pdf, BxDFType *sampledType) const {
    // Cosine-sample the hemisphere, flipping the direction if necessary
    *wi = CosineSampleHemisphere(u);
    if (wo.z < 0) wi->z *= -1;
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

Float BxDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
}

Spectrum LambertianTransmission::Sample_f(const Vector3f &wo, Vector3f *wi,
                                          const Point2f &u, Float *pdf,
                                          BxDFType *sampledType) const {
    *wi = CosineSampleHemisphere(u);
    if (wo.z > 0) wi->z *= -1;
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

Float LambertianTransmission::Pdf(const Vector3f &wo,
                                  const Vector3f &wi) const {
    return !SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
}

Spectrum MicrofacetReflection::Sample_f(const Vector3f &wo, Vector3f *wi,
                                        const Point2f &u, Float *pdf,
                                        BxDFType *sampledType) const {
    // Sample microfacet orientation $\wh$ and reflected direction $\wi$
    if (wo.z == 0) return 0.;
    Vector3f wh = distribution->Sample_wh(wo, u);
    *wi = Reflect(wo, wh);
    if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);

    // Compute PDF of _wi_ for microfacet reflection
    *pdf = distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
    return f(wo, *wi);
}

Float MicrofacetReflection::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (!SameHemisphere(wo, wi)) return 0;
    Vector3f wh = Normalize(wo + wi);
    return distribution->Pdf(wo, wh) / (4 * Dot(wo, wh));
}

Spectrum MicrofacetTransmission::Sample_f(const Vector3f &wo, Vector3f *wi,
                                          const Point2f &u, Float *pdf,
                                          BxDFType *sampledType) const {
    if (wo.z == 0) return 0.;
    Vector3f wh = distribution->Sample_wh(wo, u);
    Float eta = CosTheta(wo) > 0 ? (etaA / etaB) : (etaB / etaA);
    if (!Refract(wo, (Normal3f)wh, eta, wi)) return 0;
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

Float MicrofacetTransmission::Pdf(const Vector3f &wo,
                                  const Vector3f &wi) const {
    if (SameHemisphere(wo, wi)) return 0;
    // Compute $\wh$ from $\wo$ and $\wi$ for microfacet transmission
    Float eta = CosTheta(wo) > 0 ? (etaB / etaA) : (etaA / etaB);
    Vector3f wh = Normalize(wo + wi * eta);

    // Compute change of variables _dwh\_dwi_ for microfacet transmission
    Float sqrtDenom = Dot(wo, wh) + eta * Dot(wi, wh);
    Float dwh_dwi =
        std::abs((eta * eta * Dot(wi, wh)) / (sqrtDenom * sqrtDenom));
    return distribution->Pdf(wo, wh) * dwh_dwi;
}

Spectrum FresnelBlend::Sample_f(const Vector3f &wo, Vector3f *wi,
                                const Point2f &uOrig, Float *pdf,
                                BxDFType *sampledType) const {
    Point2f u = uOrig;
    if (u[0] < .5) {
        u[0] = std::min(2 * u[0], OneMinusEpsilon);
        // Cosine-sample the hemisphere, flipping the direction if necessary
        *wi = CosineSampleHemisphere(u);
        if (wo.z < 0) wi->z *= -1;
    } else {
        u[0] = std::min(2 * (u[0] - .5f), OneMinusEpsilon);
        // Sample microfacet orientation $\wh$ and reflected direction $\wi$
        Vector3f wh = distribution->Sample_wh(wo, u);
        *wi = Reflect(wo, wh);
        if (!SameHemisphere(wo, *wi)) return Spectrum(0.f);
    }
    *pdf = Pdf(wo, *wi);
    return f(wo, *wi);
}

Float FresnelBlend::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    if (!SameHemisphere(wo, wi)) return 0;
    Vector3f wh = Normalize(wo + wi);
    Float pdf_wh = distribution->Pdf(wo, wh);
    return .5f * (AbsCosTheta(wi) * InvPi + pdf_wh / (4 * Dot(wo, wh)));
}

Spectrum FresnelSpecular::Sample_f(const Vector3f &wo, Vector3f *wi,
                                   const Point2f &u, Float *pdf,
                                   BxDFType *sampledType) const {
    Float F = FrDielectric(CosTheta(wo), etaA, etaB);
    if (u[0] < F) {
        // Compute specular reflection for _FresnelSpecular_

        // Compute perfect specular reflection direction
        *wi = Vector3f(-wo.x, -wo.y, wo.z);
        if (sampledType)
            *sampledType = BxDFType(BSDF_SPECULAR | BSDF_REFLECTION);
        *pdf = F;
        return F * R / AbsCosTheta(*wi);
    } else {
        // Compute specular transmission for _FresnelSpecular_

        // Figure out which $\eta$ is incident and which is transmitted
        bool entering = CosTheta(wo) > 0;
        Float etaI = entering ? etaA : etaB;
        Float etaT = entering ? etaB : etaA;

        // Compute ray direction for specular transmission
        if (!Refract(wo, Faceforward(Normal3f(0, 0, 1), wo), etaI / etaT, wi))
            return 0;
        Spectrum ft = T * (1 - F);

        // Account for non-symmetry with transmission to different medium
        if (mode == TransportMode::Radiance)
            ft *= (etaI * etaI) / (etaT * etaT);
        if (sampledType)
            *sampledType = BxDFType(BSDF_SPECULAR | BSDF_TRANSMISSION);
        *pdf = 1 - F;
        return ft / AbsCosTheta(*wi);
    }
}

std::string FresnelSpecular::ToString() const {
    return std::string("[ FresnelSpecular R: ") + R.ToString() +
           std::string(" T: ") + T.ToString() +
           StringPrintf(" etaA: %f etaB: %f ", etaA, etaB) +
           std::string(" mode : ") +
           (mode == TransportMode::Radiance ? std::string("RADIANCE")
                                            : std::string("IMPORTANCE")) +
           std::string(" ]");
}

Spectrum FourierBSDF::Sample_f(const Vector3f &wo, Vector3f *wi,
                               const Point2f &u, Float *pdf,
                               BxDFType *sampledType) const {
    // Sample zenith angle component for _FourierBSDF_
    Float muO = CosTheta(wo);
    Float pdfMu;
    Float muI = SampleCatmullRom2D(bsdfTable.nMu, bsdfTable.nMu, bsdfTable.mu,
                                   bsdfTable.mu, bsdfTable.a0, bsdfTable.cdf,
                                   muO, u[1], nullptr, &pdfMu);

    // Compute Fourier coefficients $a_k$ for $(\mui, \muo)$

    // Determine offsets and weights for $\mui$ and $\muo$
    int offsetI, offsetO;
    Float weightsI[4], weightsO[4];
    if (!bsdfTable.GetWeightsAndOffset(muI, &offsetI, weightsI) ||
        !bsdfTable.GetWeightsAndOffset(muO, &offsetO, weightsO))
        return Spectrum(0.f);

    // Allocate storage to accumulate _ak_ coefficients
    Float *ak = ALLOCA(Float, bsdfTable.mMax * bsdfTable.nChannels);
    memset(ak, 0, bsdfTable.mMax * bsdfTable.nChannels * sizeof(Float));

    // Accumulate weighted sums of nearby $a_k$ coefficients
    int mMax = 0;
    for (int b = 0; b < 4; ++b) {
        for (int a = 0; a < 4; ++a) {
            // Add contribution of _(a, b)_ to $a_k$ values
            Float weight = weightsI[a] * weightsO[b];
            if (weight != 0) {
                int m;
                const Float *ap = bsdfTable.GetAk(offsetI + a, offsetO + b, &m);
                mMax = std::max(mMax, m);
                for (int c = 0; c < bsdfTable.nChannels; ++c)
                    for (int k = 0; k < m; ++k)
                        ak[c * bsdfTable.mMax + k] += weight * ap[c * m + k];
            }
        }
    }

    // Importance sample the luminance Fourier expansion
    Float phi, pdfPhi;
    Float Y = SampleFourier(ak, bsdfTable.recip, mMax, u[0], &pdfPhi, &phi);
    *pdf = std::max((Float)0, pdfPhi * pdfMu);

    // Compute the scattered direction for _FourierBSDF_
    Float sin2ThetaI = std::max((Float)0, 1 - muI * muI);
    Float norm = std::sqrt(sin2ThetaI / Sin2Theta(wo));
    if (std::isinf(norm)) norm = 0;
    Float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
    *wi = -Vector3f(norm * (cosPhi * wo.x - sinPhi * wo.y),
                    norm * (sinPhi * wo.x + cosPhi * wo.y), muI);

    // Mathematically, wi will be normalized (if wo was). However, in
    // practice, floating-point rounding error can cause some error to
    // accumulate in the computed value of wi here. This can be
    // catastrophic: if the ray intersects an object with the FourierBSDF
    // again and the wo (based on such a wi) is nearly perpendicular to the
    // surface, then the wi computed at the next intersection can end up
    // being substantially (like 4x) longer than normalized, which leads to
    // all sorts of errors, including negative spectral values. Therefore,
    // we normalize again here.
    *wi = Normalize(*wi);

    // Evaluate remaining Fourier expansions for angle $\phi$
    Float scale = muI != 0 ? (1 / std::abs(muI)) : (Float)0;
    if (mode == TransportMode::Radiance && muI * muO > 0) {
        float eta = muI > 0 ? 1 / bsdfTable.eta : bsdfTable.eta;
        scale *= eta * eta;
    }

    if (bsdfTable.nChannels == 1) return Spectrum(Y * scale);
    Float R = Fourier(ak + 1 * bsdfTable.mMax, mMax, cosPhi);
    Float B = Fourier(ak + 2 * bsdfTable.mMax, mMax, cosPhi);
    Float G = 1.39829f * Y - 0.100913f * B - 0.297375f * R;
    Float rgb[3] = {R * scale, G * scale, B * scale};
    return Spectrum::FromRGB(rgb).Clamp();
}

Float FourierBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    // Find the zenith angle cosines and azimuth difference angle
    Float muI = CosTheta(-wi), muO = CosTheta(wo);
    Float cosPhi = CosDPhi(-wi, wo);

    // Compute luminance Fourier coefficients $a_k$ for $(\mui, \muo)$
    int offsetI, offsetO;
    Float weightsI[4], weightsO[4];
    if (!bsdfTable.GetWeightsAndOffset(muI, &offsetI, weightsI) ||
        !bsdfTable.GetWeightsAndOffset(muO, &offsetO, weightsO))
        return 0;
    Float *ak = ALLOCA(Float, bsdfTable.mMax);
    memset(ak, 0, bsdfTable.mMax * sizeof(Float));
    int mMax = 0;
    for (int o = 0; o < 4; ++o) {
        for (int i = 0; i < 4; ++i) {
            Float weight = weightsI[i] * weightsO[o];
            if (weight == 0) continue;

            int order;
            const Float *coeffs =
                bsdfTable.GetAk(offsetI + i, offsetO + o, &order);
            mMax = std::max(mMax, order);

            for (int k = 0; k < order; ++k) ak[k] += *coeffs++ * weight;
        }
    }

    // Evaluate probability of sampling _wi_
    Float rho = 0;
    for (int o = 0; o < 4; ++o) {
        if (weightsO[o] == 0) continue;
        rho +=
            weightsO[o] *
            bsdfTable.cdf[(offsetO + o) * bsdfTable.nMu + bsdfTable.nMu - 1] *
            (2 * Pi);
    }
    Float Y = Fourier(ak, mMax, cosPhi);
    return (rho > 0 && Y > 0) ? (Y / rho) : 0;
}


Spectrum BxDF::rho(const Vector3f &w, int nSamples, const Point2f *u) const {
    Spectrum r(0.);
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hd}$
        Vector3f wi;
        Float pdf = 0;
        Spectrum f = Sample_f(w, &wi, u[i], &pdf);
        if (pdf > 0) r += f * AbsCosTheta(wi) / pdf;
    }
    return r / nSamples;
}

Spectrum BxDF::rho(int nSamples, const Point2f *u1, const Point2f *u2) const {
    Spectrum r(0.f);
    for (int i = 0; i < nSamples; ++i) {
        // Estimate one term of $\rho_\roman{hh}$
        Vector3f wo, wi;
        wo = UniformSampleHemisphere(u1[i]);
        Float pdfo = UniformHemispherePdf(), pdfi = 0;
        Spectrum f = Sample_f(wo, &wi, u2[i], &pdfi);
        if (pdfi > 0)
            r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdfo * pdfi);
    }
    return r / (Pi * nSamples);
}

// BSDF Method Definitions
Spectrum BSDF::f(const Vector3f &woW, const Vector3f &wiW,
                 BxDFType flags) const {
    ProfilePhase pp(Prof::BSDFEvaluation);
    Vector3f wi = WorldToLocal(wiW), wo = WorldToLocal(woW);
    if (wo.z == 0) return 0.;
    /*
    if (fabs(wo.z) > 1 || fabs(wi.z) > 1) {
        std::cout<< " wo " << wo << "\n";
        std::cout<< " wi " << wi << "\n";
        std::cout<< " ss " << ss << "\n";
        std::cout<< " ts " << ts << "\n";
        std::cout<< " ns " << ns << "\n";
    }
    */
    bool reflect = Dot(wiW, ng) * Dot(woW, ng) > 0;
    Spectrum f(0.f);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags) &&
            ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
             (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
            f += bxdfs[i]->f(wo, wi);
    return f;
}

Spectrum BSDF::rho(int nSamples, const Point2f *samples1,
                   const Point2f *samples2, BxDFType flags) const {
    Spectrum ret(0.f);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            ret += bxdfs[i]->rho(nSamples, samples1, samples2);
    return ret;
}

Spectrum BSDF::rho(const Vector3f &wo, int nSamples, const Point2f *samples,
                   BxDFType flags) const {
    Spectrum ret(0.f);
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags))
            ret += bxdfs[i]->rho(wo, nSamples, samples);
    return ret;
}

Spectrum BSDF::Sample_f(const Vector3f &woWorld, Vector3f *wiWorld,
                        const Point2f &u, Float *pdf, BxDFType type,
                        BxDFType *sampledType) const {
    ProfilePhase pp(Prof::BSDFSampling);
    // Choose which _BxDF_ to sample
    int matchingComps = NumComponents(type);
    if (matchingComps == 0) {
        *pdf = 0;
        if (sampledType) *sampledType = BxDFType(0);
        return Spectrum(0);
    }
    int comp =
        std::min((int)std::floor(u[0] * matchingComps), matchingComps - 1);

    // Get _BxDF_ pointer for chosen component
    BxDF *bxdf = nullptr;
    int count = comp;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(type) && count-- == 0) {
            bxdf = bxdfs[i];
            break;
        }
    CHECK_NOTNULL(bxdf);
    VLOG(2) << "BSDF::Sample_f chose comp = " << comp << " / matching = " <<
        matchingComps << ", bxdf: " << bxdf->ToString();

    // Remap _BxDF_ sample _u_ to $[0,1)^2$
    Point2f uRemapped(std::min(u[0] * matchingComps - comp, OneMinusEpsilon),
                      u[1]);

    // Sample chosen _BxDF_
    Vector3f wi, wo = WorldToLocal(woWorld);
    if (wo.z == 0) return 0.;
    *pdf = 0;
    if (sampledType) *sampledType = bxdf->type;
    Spectrum f = bxdf->Sample_f(wo, &wi, uRemapped, pdf, sampledType);
    VLOG(2) << "For wo = " << wo << ", sampled f = " << f << ", pdf = "
            << *pdf << ", ratio = " << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.))
            << ", wi = " << wi;
    if (*pdf == 0) {
        if (sampledType) *sampledType = BxDFType(0);
        return 0;
    }
    *wiWorld = LocalToWorld(wi);

    // Compute overall PDF with all matching _BxDF_s
    if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1)
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i] != bxdf && bxdfs[i]->MatchesFlags(type))
                *pdf += bxdfs[i]->Pdf(wo, wi);
    if (matchingComps > 1) *pdf /= matchingComps;

    // Compute value of BSDF for sampled direction
    if (!(bxdf->type & BSDF_SPECULAR) && matchingComps > 1) {
        bool reflect = Dot(*wiWorld, ng) * Dot(woWorld, ng) > 0;
        f = 0.;
        for (int i = 0; i < nBxDFs; ++i)
            if (bxdfs[i]->MatchesFlags(type) &&
                ((reflect && (bxdfs[i]->type & BSDF_REFLECTION)) ||
                 (!reflect && (bxdfs[i]->type & BSDF_TRANSMISSION))))
                f += bxdfs[i]->f(wo, wi);
    }
    VLOG(2) << "Overall f = " << f << ", pdf = " << *pdf << ", ratio = "
            << ((*pdf > 0) ? (f / *pdf) : Spectrum(0.));
    return f;
}

Float BSDF::Pdf(const Vector3f &woWorld, const Vector3f &wiWorld,
                BxDFType flags) const {
    ProfilePhase pp(Prof::BSDFPdf);
    if (nBxDFs == 0.f) return 0.f;
    Vector3f wo = WorldToLocal(woWorld), wi = WorldToLocal(wiWorld);
    if (wo.z == 0) return 0.;
    Float pdf = 0.f;
    int matchingComps = 0;
    for (int i = 0; i < nBxDFs; ++i)
        if (bxdfs[i]->MatchesFlags(flags)) {
            ++matchingComps;
            pdf += bxdfs[i]->Pdf(wo, wi);
        }
    Float v = matchingComps > 0 ? pdf / matchingComps : 0.f;
    return v;
}

std::string BSDF::ToString() const {
    std::string s = StringPrintf("[ BSDF eta: %f nBxDFs: %d", eta, nBxDFs);
    for (int i = 0; i < nBxDFs; ++i)
        s += StringPrintf("\n  bxdfs[%d]: ", i) + bxdfs[i]->ToString();
    return s + std::string(" ]");
}

//////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
////GaussianBSDF:  Mandy and Feng ////////////////////////////////////////////////////////////
////3D fitting of both single and multiple scattering/////////////////////////////////////////
////Has artifacts don't use //////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

//GaussianBSDF methods
//Uniform sampling for now for Sample_f and Pdf computation
//TODO: add importance sampling of GMM
Spectrum GaussianBSDF::Sample_f(const Vector3f &wo, Vector3f *wi,
                                          const Point2f &u, Float *pdf,
                                          BxDFType *sampledType) const {
    //return MicrofacetReflection::Sample_f(wo, wi, u, pdf, sampledType);

    //Uniform sampling for gaussian bsdf 
    Vector3f wo_normalized = Normalize(wo);

    // Sample gaussian orientation $\wh$ and reflected direction $\wi$
    *wi = CosineSampleHemisphere(u);
    *wi = Normalize(*wi);
    if (wo.z < 0) wi->z *= -1;
    if (!SameHemisphere(wo_normalized, *wi) || wi->z == 0) return Spectrum(0.f);
    *pdf = Pdf(wo_normalized, *wi);
    return f(wo_normalized, *wi);
}

Float GaussianBSDF::Pdf(const Vector3f &wo, const Vector3f &wi) const {
    //return MicrofacetReflection::Pdf(wo, wi);
    
    //Uniform sampling for gaussian bsdf 
    return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * InvPi : 0;
}

/*
//uncomment to debug for negative wo and wi
void check_omega_oi_negative(const Vector3f &wo, const Vector3f &wi) {
    std::cout<< "wo " << wo<< " wi " << wi << "\n";
    fflush(stdout);
}
*/

Spectrum GaussianBSDF::f(const Vector3f &woO, const Vector3f &wiO) const {

    //This evaluation is here only for potential debugging    
    //Spectrum mf =  MicrofacetReflection::f(woO, wiO);

    Vector3f wo = Normalize(woO);
    Vector3f wi = Normalize(wiO);

    if (wo.z * wi.z < 0) { 
        return Spectrum(0.);
    }

    if (wo.z < 1e-6 || wi.z < 1e-6) {
        return 0;
    }
    Float cosThetaI = AbsCosTheta(wi);
    Float cosThetaO = AbsCosTheta(wo);
    // handle degenerate cases
    if (cosThetaI == 0 || cosThetaO == 0) return Spectrum(0.);

    Float phi = atan2(wo.y, wo.x);
    if (phi < 0) phi += 2.0 * M_PI;
    if (fabs(phi) > 1e-5) {
        phi = Degrees(phi);
        Transform rotation = RotateZ(-phi);
        wo = rotation(wo);
        if (fabs(wo.y) > 1e-3) {
            std::cout << "phi: " << phi << "mu: " << cosThetaO << "\n";
            std::cout<< "woO " << woO<< " wo " << wo << "\n";
            std::cout<< "rotation" << rotation << "\n";
            fflush(stdout);
        }
        assert(fabs(wo.y) < 1e-5);
        wi = rotation(wi);
    }

    wo.z = abs(wo.z);
    wi.z = abs(wi.z);
   
    Float J = computeJacobian(wi, wo, 1, 1);
    
    Vector3f wh = wi + wo;
    wh = Normalize(wh);
    
    // calculate probability in slope domain using fitted GMM
    Float x = wh.x/wh.z;
    Float y = wh.y/wh.z;
    Float zo = wo.z > 1? 1: wo.z;
    zo = acos(zo); 
    if (zo <0 || zo >= M_PI * .5 || isNaN(zo)) {
        std::cout<< "woO " << woO<< " wo " << wo << "\n";
        std::cout << "gaussian brdf value: 0" << "\n";
        fflush(stdout);
        return Spectrum(0);
    }
    Float p = gm->prob(x,y,zo);
    if (isNaN(p)) {
        std::cout<< "has NaN prob" << "\n";
        fflush(stdout);
        return Spectrum(0);
    }
    Spectrum brdf =  p * J / cosThetaI;

    /*
    if (mf[1] > 2.0 * brdf[1]) {
        std::cout << "beckmann brdf value:" << mf << "\n";
        std::cout << "gaussian brdf value:" << brdf << "\n";
        fflush(stdout);
    }
    */
    return brdf;
}

std::string GaussianBSDF::ToString() const {
     return MicrofacetReflection::ToString();
     //TODO: add gaussian mixture string rep
}



}  // namespace pbrt
