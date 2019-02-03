#include "gaussianscatter.h"

/*
    -0.0433   -0.0395    0.0097    0.2321
    0.0552   -0.0308    0.2410    0.0070
    0.0168   -0.0157    0.0013   -0.0024
    0.0194   -0.0319    0.0027    0.0626
   -0.0100    0.0759   -0.0081    0.0631
*/

using namespace pbrt;

GaussianScatter::GaussianScatter() : energyOnly(false) {
    penergy = Polynomial(-0.0433,   -0.0395,    0.0097,    0.2321);
    pmean[0] = Polynomial(0.0552,   -0.0308,    0.2410,   0.0070);
    pmean[1] = Polynomial( 0.0168,   -0.0157,    0.0013,   -0.0024);
    pcov[0] = Polynomial(0.0194,   -0.0319,    0.0027,    0.0626);
    pcov[1] = Polynomial(-0.0100,    0.0759,   -0.0081,    0.0631);

    validScatter = true;
    alpha = .7;
    validateScatter();
}

GaussianScatter::GaussianScatter(Float alpha, bool energyOnly, bool isTransmission): alpha(alpha), 
                 energyOnly(energyOnly), isTransmission(isTransmission) {

    std::cout << "energy only: " << energyOnly <<"\n";
    char buf[10];
    sprintf(buf, "%1.1f", alpha);
    std::string prefix_dir = "MultiScatteringCoef/";
    std::string prefix = isTransmission? "transmission_coef/" : "reflection_coef/";
    std::string scatterName = prefix_dir + prefix + 
            "scatter_" + std::string(buf) + ".txt";
    std::cout << "load scatter file:" << scatterName << "\n";
    std::ifstream scatterfile(scatterName);
    validScatter = false;

    if (scatterfile.is_open()) {
        int coefCount = 0;
        scatterfile >> coefCount;
        std::cout<<" coef count for scatter "<< coefCount << "\n";
        validScatter = coefCount > 0;
        validScatter &= penergy.readFromFile(scatterfile, coefCount);
        validScatter &= pmean[0].readFromFile(scatterfile, coefCount);
        validScatter &= pmean[1].readFromFile(scatterfile, coefCount);
        validScatter &= pcov[0].readFromFile(scatterfile, coefCount);
        validScatter &= pcov[1].readFromFile(scatterfile, coefCount);
        validateScatter();
    }
}

bool Polynomial::readFromFile(std::ifstream& file, int coefCount) {

    coef.resize(coefCount);
    for (int i = 0; i < coefCount; i++) {
        file >> coef[i];
        std::cout << coef[i];
    }
    std::cout << "\n";
    
    return file.good();
}

void GaussianScatter::validateScatter() {
    if (validScatter) {
        Polynomial emu = penergy.multiply_x();
        Polynomial ienergy = emu.integral();
        energyAve = ienergy.eval(1.0) - ienergy.eval(0);
        energyAve = 2.0 * energyAve;
        
        Float energy = penergy.eval(.5);
        Float m0 = pmean[0].eval(.5);
        Float m1 = pmean[1].eval(.5);
        Float c0 = pcov[0].eval(.5);
        Float c1 = pcov[1].eval(.5);
        std::cout << "energy, m0, m1, c0, c1, energyAve "<< energy<< " "<< m0<< " "<< m1<< " "<< c0<< " "<< c1 << " " << energyAve <<"\n";
    }
}

