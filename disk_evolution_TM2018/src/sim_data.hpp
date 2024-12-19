#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP

#include "array.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "define_macro.hpp"
#include "grid.hpp"

#include <iostream>
#include <cmath>



constexpr double SIGMA_GAS_FLOOR  = 1.0e-15;
constexpr double SIGMA_DUST_FLOOR = 1.0e-15;

// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

class Gas
{
public:
    Gas();
    ~Gas();

    void SetGas(Utils::InputConfigure &input, Grid &grid);

    array::Double1D Sigma_;          // surface density [g cm^-2]
    array::Double1D Sigma_floor_;    // floor value of surface density [g cm^-2]
    array::Double1D T_;              // temperature [K]
    array::Double1D cs_;             // sound speed [cm s^-1]
    array::Double1D Omega_;          // angular velocity [s^-1]
    array::Double1D j_;              // specific angular moemntum [cm^2 s^-1]
    array::Double1D Hg_;             // gas scale height [cm]
    array::Double1D rho_mid_;        // gas mass density [g cm^-3]
    array::Double1D P_;              // gas pressure at disk midplane []
    array::Double1D P_integ_;        // vertically integrated gas pressure []
    array::Double1D Qt_;             // Toomre's Q parameter
    array::Double1D Mr_bnd_;         // enclosed mass at cell boundary [g]
    array::Double1D Mr_cen_;         // enclosed mass at cell center [g]
    array::Double1D alpha_;          // alpha parameter
    array::Double1D vr_;             // radial velocity [cm s^-1]
    array::Double1D vr_vis_;         // viscous velocity [cm s^-1]
    array::Double1D vr_eta_;         // gas pressure gradient velocity [cm s^-1]
    array::Double1D vr_src_;         // source term velocity [cm s^-1]
    array::Double2D cvg_;            // gas velocity coefficient by gas and dust interaction
    array::Double1D Sigma_dot_inf_;  // [g cm^-2 s^-1]
    array::Double1D Sigma_dot_wind_; // [g cm^-2 s^-1]
    double alpha_turb_;              // alpha parameter of turbulence (MRI)
    double c_wind_;                  // wind strength
    double total_mass_;              // total disk gas mass [g]
private:
    bool is_allocate_arrays_;
};


// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

class Dust
{
public:
    Dust();
    ~Dust();

    void SetDust(Utils::InputConfigure &input, Grid &grid);

    array::Double1D Sigma_;         // dust surface density
    array::Double1D Sigma_floor_;   // floor value of dust surface density
    array::Double1D vr_;            // dust radial velocity
    array::Double2D cvd_;           // dust velocity coefficient by gas and dust interaction
    double fdg_;                    // dust to gas mass ratio in envelope
    double St_;                     // Stokes number 
    double total_mass_;             // total disk dust mass
private:
    bool is_allocate_arrays_;
};

// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

class Cloud
{
public:
    void SetCloud(Utils::InputConfigure &input);
    int nshells_;      // number of cloud shells
    double nc_;        // central number density
    double T_;         // temperature
    double Omega_c_;   // angular velocity
    double f_;         // gravity enhancement factor
    double Mc_;        // total cloud mass
    double e_env_;     // internal energy of envelope
    double rini_;      // initial radius of infall material
    double drinidt_;   // 
    double mdot_inf_;  // mass infall rate from envelope to star-disk system
private:
};

// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

class Star
{
public:
    void SetStar(Utils::InputConfigure &input);
    double mass_init_;    // initial star mass
    double mass_;         // star mass
    double Ls_;           // star luminosity
    double Lacc_;         // accretion luminosity
private:
};

// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

class Time
{
public:
    void SetTime(Utils::InputConfigure &input);
    double tend_;
    double t_;
    double cfl_;
    double delta_tout_;
private:
};

// ################################################################################## //
// ################################################################################## //
// ################################################################################## //

class SimulationData
{
public:

    SimulationData(Utils::InputConfigure &input);
    ~SimulationData();

    Grid grid_;
    Gas  gas_;
    Dust dust_;
    Cloud cloud_;
    Star star_;
    Time time_;

    std::string outdir_;
    int output_count_;

    double total_disk_mass_;
    double total_mass_;             // star + disk mass
    double mdot_acc_disk_;          // 内側境界から流出する量 (mass accretion rate from disk to star)
    double mdot_acc_env_;           // mass accretion rate from envelope to satr
    double mdot_inf_;               // infall rate from envelope to disk
    double total_infall_mass_;
    double mdot_wind_;              // wind mass loss rate from disk
    double mdot_wind_sweep_;
    double wind_loss_mass_;
    double wind_loss_mass_sweep_;
    double total_wind_loss_mass_;

    void OutputGrid();
    void OutputData(const std::string filename);
    void LogOutput(const int count, std::ofstream &file);

private:

    void Initialization();

};

#endif /* SIM_DATA_HPP */