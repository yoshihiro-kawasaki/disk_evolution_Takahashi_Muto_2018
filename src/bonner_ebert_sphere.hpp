#ifndef _BONNER_EBERT_SPHERE_HPP_
#define _BONNER_EBERT_SPHERE_HPP_

#include "array.hpp"
#include "constants.hpp"
#include "define_macro.hpp"
#include "sim_data.hpp"

class BonnerEbertSphere
{
public:

    BonnerEbertSphere(SimulationData *pdata_);
    ~BonnerEbertSphere();

    void SetBonnerEbertSphere();
    inline double GetRoutBonner() { return rout_bonner_; };
    double CalculateInitialRinit();
    void CalculateMdotInfall(const double r_init, double &drinitdt, double &mdot_inf);
    void Output(const std::string file_name);

protected:

    SimulationData *pdata_;

    int    nshells_;                // number of cloud shells
    double Tc_;                    // BE-sphere (cloud) temperature [K]
    double cs_c_;                  // sound speed of BE-sphere [cm s^-1]
    double Mc_;                    // total BE-sphere mass [g]
    double nc_;                    // central number density of BE-sphere [cm^-3] 
    double rho_c_;                 // central mass density of BE-sphere [g cm^-3]
    double Omega_c_;               // angular speed of BE-sphere [s^-1]
    double f_;                     // density factor to promote gravitational collapse takahashi et al 2013
    double Mp_init_;               // initial protostellar mass [g]
    double rout_bonner_;           // outer edge of BE-sphere (cloud) [cm]
    double rn_out_bonner_;         // nomalized outer edge of BE-sphere (cloud)
    
    array::Double1D rn_;           // normalized length
    array::Double1D psi_;          // normalized gravitatio potential
    array::Double1D dpsidrn_;      // derivative of psi_
    array::Double1D rho_bonner_;   // cloud mass density [g/cm^3]
    array::Double1D r_bonner_;     // cloud radius [cm]
    array::Double1D mass_bonner_;  // cloud mass [g]
    
    int ib_con_;                   // current index of shells

private:

    void SolveLaneEmdenEquation(const int n, double &x, const double h, array::Double1D y);

};

#endif /* _BONNER_EBERT_SPHERE_HPP_ */