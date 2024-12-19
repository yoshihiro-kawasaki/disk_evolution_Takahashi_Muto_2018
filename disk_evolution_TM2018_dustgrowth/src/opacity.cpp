#include "opacity.hpp"

double OpacityModelZhu2012(const double P, const double T)
{
    // double P = rhog * T * pc::BOLTZMANN_CONSTANT / pc::GAS_MOLECULAR_MASS;
    double logP = std::log10(P);
    double logT = std::log10(T);
    double logkappa;

    if (logT < 0.01899*logP + 2.14072) {
        // Water ice opacity
        logkappa = 1.5*logT - 2.481329;
    } 
    else if (logT < 0.01899*logP + 2.21218) {
        // Water ice evaporation
        logkappa = -3.53154212*logT + 0.094991625*logP + 8.292767875; 
    } 
    else if (logT < 2.79055) {
        // Metal grain opacity
          logkappa = 1.5*logT -2.840695;
    } 
    else if (logT < 2.96931) {
        // Graphite corrosion
        logkappa = -5.832*logT + 17.7;
    } 
    else if (logT < 0.0268175*logP + 3.1569625) {
        // Grain opacity
        logkappa = 2.129*logT - 5.9398;
    } 
    else if (logT < 0.028084*logP + 3.192336) {
        // Silicate grain evaporation
        logkappa = -42.98075*logT + 1.3115765*logP + 135.127016;
    } 
    else if (logT < 0.03*logP + 3.28) {
        // Water vapor
        logkappa = 4.0625*logT - 15.0125;
    } 
    else if (logT < 0.00832*logP + 3.41) {
        logkappa = -18.4808*logT + 0.6763*logP + 58.9294;
    } 
    else if (logT < 0.015*logP + 3.7) {
        // Molecular opacities
        logkappa = 2.90477*logT + 0.498325*logP - 13.9953;
    } 
    else if (logT < 0.04*logP + 3.91) {
        // H scattering
        logkappa = 10.1935*logT + 0.3821*logP -40.9361;
    } 
    else if (logT < 0.2797*logP + 3.6933) {
        // Bound-free, free-free
        logkappa = -3.3647*logT + 0.92795*logP + 12.0258;
    } 
    else {
        // Electron scattering
        logkappa = -0.48;
    }


    if ((logkappa < 3.586*logT - 16.85) && (logT < 4.0)) {
        logkappa = 3.586*logT - 16.85;
    }

    return std::pow(10.0, logkappa);
}


double OpacityModelSuzukiPlusZhu(const double P, const double T)
{
    // double P = rhog * T * pc::BOLTZMANN_CONSTANT / pc::GAS_MOLECULAR_MASS;
    double logP = std::log10(P);
    double logT = std::log10(T);
    double kappad, kappag, logkappa;

    kappad = 2.25 * (1.0 - std::tanh((T -1.5e3)*0.006666666666666667)) * std::min(1.0, SQR(T *0.006666666666666667));
    kappag = 0.0;

    if (logT < 0.03*logP + 3.28) {
        // Water vapor
        logkappa = 4.0625*logT - 15.0125;
    } 
    else if (logT < 0.00832*logP + 3.41) {
        logkappa = -18.4808*logT + 0.6763*logP + 58.9294;
    } 
    else if (logT < 0.015*logP + 3.7) {
        // Molecular opacities
        logkappa = 2.90477*logT + 0.498325*logP - 13.9953;
    } 
    else if (logT < 0.04*logP + 3.91) {
        // H scattering
        logkappa = 10.1935*logT + 0.3821*logP -40.9361;
    } 
    else if (logT < 0.2797*logP + 3.6933) {
        // Bound-free, free-free
        logkappa = -3.3647*logT + 0.92795*logP + 12.0258;
    } 
    else {
        // Electron scattering
        logkappa = -0.48;
    }


    if ((logkappa < 3.586*logT - 16.85) && (logT < 4.0)) {
        logkappa = 3.586*logT - 16.85;
    }

    kappag = std::pow(10.0, logkappa);

    return kappad + kappag;
}

