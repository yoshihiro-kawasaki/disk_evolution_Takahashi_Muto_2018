#ifndef SIM_DATA_HPP
#define SIM_DATA_HPP

#include "array.hpp"
#include "utils.hpp"
#include "physical_constants.hpp"
#include "define_macro.hpp"


// constant
constexpr double sqrt_2pi = 2.5066282746310002; // sqrt(2.0*pi)

constexpr double inv_AU = 1.0 / pc::ASTRONOMICAL_UNIT; // 1.0 / 1au

constexpr double inv_mg = 1.0 / pc::GAS_MOLECULAR_MASS; // 1.0 / mg

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
        Array1D<double> dr_;
        Array1D<double> inv_r_cen_;
        Array1D<double> inv_r_vol_;
        std::string spacing_;
    } grid_;

    struct Gas
    {
        Array1D<double> Sigma_;        // surface density
        Array1D<double> Sigma_floor_;
        Array1D<double> T_;            // temperature
        Array1D<double> cs_;           // sound speed
        Array1D<double> Omega_;        // angular velocity
        Array1D<double> j_;
        Array1D<double> Hg_;
        Array1D<double> P_;
        Array1D<double> Mr_;           // enclosed mass
        Array1D<double> alpha_;        // alpha parameter
        Array1D<double> vr_;
        Array1D<double> vr_vis_;
        Array1D<double> vr_eta_;
        Array1D<double> vr_src_;
        double alpha_turb_;
        double c_wind_;
        double total_mass_;
    } gas_;

    struct Dust
    {
        Array1D<double> Sigma_;
        Array1D<double> vr_;
        double fdg_;
        double St_;
        double total_mass_;
    } dust_;

    struct Cloud
    {
        int nshells_;      // number of cloud shells
        double nc_;        // central number density
        double T_;         // temperature
        double Omega_c_;   // angular velocity
        double f_;         // 
        double Mc_;        // total cloud mass
    } cloud_;


    struct Star
    {
        double mass_init_;    // initial star mass
        double mass_;        // star mass
    } star_;

    double time_;

    double mdot_;
    double total_disk_mass_;
    double total_mass_;
    double total_wind_loss_mass_;

    std::string outdir_;

    void Output(std::string filename);

private:

    void SetLogSpacingGrid();
    void SetRootSpacingGrid();

};

#endif /* SIM_DATA_HPP */