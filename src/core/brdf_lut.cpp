#include <fstream>
#include <sstream>
#include "brdf_lut.h"

namespace pbrt{

std::string makeFileName(const std::string& prefix, 
      Float alpha, uint rayDepth) {

    std::ostringstream oss;
    oss << prefix << "_d_"<< rayDepth << "_a_" << alpha<<".p";
    std::cout<<"brdf filename: " << oss.str() <<"\n";
    return oss.str();
}

BrdfLUT::BrdfLUT(const std::string& filePrefix, Float alpha, uint rayDepth)            :alpha(alpha), rayDepth(rayDepth) {
    auto filename = makeFileName(filePrefix, alpha, rayDepth);
    load(filename);
}

void BrdfLUT::save(const std::string& filePrefix) {
    if (!normed) normalize();
    auto filename = makeFileName(filePrefix, alpha, rayDepth);

    std::ofstream output;
    output.open(filename, std::ios::out | std::ios::binary);
  
    output.write((char*) &energyRatio, sizeof(energyRatio)); 
    output.write((char*) &totalSamples, sizeof(totalSamples)); 
    output.write((char*) &rayDepth, sizeof(rayDepth)); 
    output.write((char*) &dim, sizeof(dim)); 
    output.write((char*) pdf, sizeof(Float)* tsize);
    output.write((char*) fresnel, sizeof(Float)* tsize * 3);
    output.close();
}

void BrdfLUT::load(const std::string& filename) {
    
    std::ifstream output;
    output.open(filename, std::ios::in | std::ios::binary);
  
    output.read((char*) &energyRatio, sizeof(energyRatio)); 
    output.read((char*) &totalSamples, sizeof(totalSamples)); 
    output.read((char*) &rayDepth, sizeof(rayDepth)); 
    output.read((char*) &dim, sizeof(dim)); 
    std::cout << "dim "<< dim << "\n";
    tsize = dim * dim * dim;
    std::cout << "tsize "<< tsize << "\n";
    pdf = new Float [tsize];
    output.read((char*) pdf, sizeof(Float)* tsize);
    std::cout << "read pdf\n";
    int pdfvalid = 0;
    for (int i = 0; i < tsize; i++) {
        if (pdf[i] > 0) pdfvalid++;
    }
    std::cout << "valid pdfs: " << pdfvalid << "\n";

    fresnel = new Float[tsize * 3];
    output.read((char*) fresnel, sizeof(Float)* tsize * 3);
    output.close();
    computeUnits();
    int fresnelvalid = 0;
    for (int i = 0; i < tsize * 3; i++) {
        if (fresnel[i] > 0) fresnelvalid++;
    }
    std::cout << "valid fresnels: " << fresnelvalid << "\n";
}

} //end namespace


