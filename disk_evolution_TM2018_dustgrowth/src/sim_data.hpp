#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP

#include "array.hpp"
#include "utils.hpp"
#include "physical_constants.hpp"
#include "define_macro.hpp"

#include <iostream>
#include <cmath>


// constant
constexpr double sqrt_2pi = 2.5066282746310002; // sqrt(2.0*pi)

constexpr double inv_AU = 1.0 / pc::ASTRONOMICAL_UNIT; // 1.0 / 1au

constexpr double inv_mg = 1.0 / pc::GAS_MOLECULAR_MASS; // 1.0 / mg

constexpr double gamma_h = 1.4;

constexpr double inv_gamma_h_1 = 1.0 / (gamma_h - 1.0);

constexpr double coe_EtoT = (gamma_h - 1.0)*pc::GAS_MOLECULAR_MASS / pc::BOLTZMANN_CONSTANT;

constexpr double inv_4pi = 1.0 / (4.0*M_PI);

constexpr double inv_3pi = 1.0 / (3.0*M_PI);

//

class SimulationData
{
public:

    SimulationData(Utils::InputConfigure &input);
    ~SimulationData();

    void Initialization();


    struct Grid 
    {
        int nr_;
        int nrtotc_;
        int nrtotb_;
        int is_;
        int ie_;
        double rmin_;
        double rmax_;
        Array1D<double> r_cen_;
        Array1D<double> r_bnd_;
        Array1D<double> r_vol_;
        Array1D<double> dr_cen_;
        Array1D<double> dr_bnd_;
        Array1D<double> inv_r_cen_;
        Array1D<double> inv_r_vol_;
        Array1D<double> inv_dr_cen_;
        Array2D<double> cdiff_;      // coefficient for numerical derivative
        std::string spacing_;
    } grid_;

    struct Gas
    {
        Array1D<double> Sigma_;          // surface density
        Array1D<double> Sigma_floor_;    // floor value of surface density
        Array1D<double> T_;              // temperature
        Array1D<double> cs_;             // sound speed
        Array1D<double> Omega_;          // angular velocity
        Array1D<double> j_;              // specific angular moemntum
        Array1D<double> Hg_;             // gas scale height
        Array1D<double> P_;              // gas pressure at disk midplane
        Array1D<double> P_integ_;        // vertically integrated gas pressure
        Array1D<double> Qt_;             // Toomre's Q parameter
        Array1D<double> Mr_;             // enclosed mass
        Array1D<double> alpha_;          // alpha parameter
        Array1D<double> vr_;             // radial velocity
        Array1D<double> vr_vis_;         // viscous velocity
        Array1D<double> vr_eta_;         // gas pressure gradient velocity
        Array1D<double> vr_src_;         // source term velocity
        Array2D<double> cvg_;            // gas velocity coefficient by gas and dust interaction
        Array1D<double> Sigma_dot_inf_;  //
        Array1D<double> Sigma_dot_wind_; //
        double alpha_turb_;
        double c_wind_;                  // wind strength
        double total_mass_;              // total disk gas mass
    } gas_;

    struct Dust
    {
        Array1D<double> Sigma_;         // dust surface density
        Array1D<double> Sigma_floor_;   // floor value of dust surface density
        Array1D<double> m_;
        Array1D<double> a_;
        Array1D<double> mSigma_;
        Array1D<double> St_;
        Array1D<double> Hd_;
        Array1D<double> vr_;            // dust radial velocity
        Array2D<double> cvd_;           // dust velocity coefficient by gas and dust interaction
        double vfrag_;
        double fdg_;                    // dust to gas mass ratio in envelope
        double rho_int_;
        double a0_;
        double m0_;
        double total_mass_;             // total disk dust mass
    } dust_;
    

    struct Cloud
    {
        int nshells_;      // number of cloud shells
        double nc_;        // central number density
        double T_;         // temperature
        double Omega_c_;   // angular velocity
        double f_;         // gravity enhancement factor
        double Mc_;        // total cloud mass
        double e_env_;     // internal energy of envelope
    } cloud_;


    struct Star
    {
        double mass_init_;    // initial star mass
        double mass_;         // star mass
        double Ls_;           // star luminosity
        double Lacc_;         // accretion luminosity
    } star_;

    double time_;

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

    std::string outdir_;

    void OutputGrid();
    void OutputData(std::string filename);
    void LogOutput(const int count, std::ofstream &file);

private:

    void SetLogSpacingGrid();
    void SetRootSpacingGrid();
    void SetNumericalDeviative();

};

#endif /* SIM_DATA_HPP */