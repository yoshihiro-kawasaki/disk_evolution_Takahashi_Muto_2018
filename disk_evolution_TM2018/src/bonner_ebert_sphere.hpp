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

    int number_of_cloud_shells_;
    int ibmax_;

    double cloud_temperature_;                // [K]
    double cloud_sound_speed_;                // [cm s^-1]
    double total_cloud_mass_;                 // [g]
    double central_number_density_of_cloud_;  // [cm^-3] 
    double central_mass_density_of_cloud_;    // [g cm^-3]
    double cloud_angular_velocity_;           // [s^-1]
    double density_factor_;                   // density factor to promote gravitational collapse takahashi et al 2013
    double initial_protostellar_mass_;        // [g]
    double rout_bonner_;                      // outer edge of cloud [cm]
    double rn_out_bonner_;                    // nomalized outer edge of cloud
    int ib_con_;
    
    array::Double1D rn_;                      // normalized length
    array::Double1D psi_;                     // normalized gravitatio potential
    array::Double1D dpsidrn_;
    array::Double1D rho_bonner_;              // cloud mass density [g/cm^3]
    array::Double1D r_bonner_;                // cloud radius [cm]
    array::Double1D mass_bonner_;             // cloud mass [g]

private:

    void SolveLaneEmdenEquation(const int n, double &x, const double h, array::Double1D y);

};

#endif /* _BONNER_EBERT_SPHERE_HPP_ */