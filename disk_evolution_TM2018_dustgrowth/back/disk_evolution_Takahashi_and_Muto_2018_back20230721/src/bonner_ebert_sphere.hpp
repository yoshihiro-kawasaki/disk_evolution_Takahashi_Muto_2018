#ifndef _BONNER_EBERT_SPHERE_HPP_
#define _BONNER_EBERT_SPHERE_HPP_

#include "array.hpp"
#include "physical_constants.hpp"
#include "define_macro.hpp"
#include "explicit_solver/explicit_solver.hpp"
#include "sim_data.hpp"

class BonnerEbertSphere
{
public:

    BonnerEbertSphere(SimulationData *pdata_);
    ~BonnerEbertSphere();

    void SetBonnerEbertSphere();

    // inline double GetRadiusBonner(const int i) { return r_bonner_(i); };
    // inline double GetRhoBonner(const int i) { return rho_bonner_(i); };
    inline double GetRoutBonner() { return rout_bonner_; };

    double CalculateInitialRinit();

    void CalculateMdotInfall(const double r_init, double &drinitdt, double &mdot_inf);

    void Output(std::string file_name);

protected:

    SimulationData *pdata_;

    int number_of_cloud_shells_;
    int ibmax_;

    double cloud_temperature_;
    double cloud_sound_speed_;
    double total_cloud_mass_;
    double central_number_density_of_cloud_;
    double central_mass_density_of_cloud_;
    double cloud_angular_velocity_;
    double density_factor_;                      // density factor to promote gravitational collapse takahashi et al 2013
    double initial_protostellar_mass_;
    double rout_bonner_;           // outer edge of cloud [cm]
    double rn_out_bonner_;         // nomalized outer edge of cloud
    int ib_con_;
    
    Array1D<double> rn_;              // normalized length
    Array1D<double> psi_;             // normalized gravitatio potential
    Array1D<double> dpsidrn_;
    Array1D<double> rho_bonner_;      // cloud mass density [g/cm^3]
    Array1D<double> r_bonner_;        // cloud radius [cm]
    Array1D<double> mass_bonner_;     // cloud mass [g]

    static void Func(double x, Array1D<double> &y, Array1D<double> &dydx, void *data);

private:

};

#endif /* _BONNER_EBERT_SPHERE_HPP_ */