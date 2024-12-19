#ifndef DISK_EVOLUTION_HPP_
#define DISK_EVOLUTION_HPP_

#include "array.hpp"
#include "define_macro.hpp"
#include "sim_data.hpp"
#include "opacity.hpp"
#include "reconstruction.hpp"

constexpr double ONE_THIRD = 1.0 / 3.0;
constexpr double SQRT2 = M_SQRT2;
constexpr double SQRT_PI = 1.7724538509055159;

constexpr double COE_SPHERE = 4.0 * M_PI / 3.0;
constexpr double COE_INV_SPHERE = 1.0 / COE_SPHERE;
constexpr double COE_VBR = 8.0 * pc::BOLTZMANN_CONSTANT / M_PI;

class DiskEvolution
{
public:

    DiskEvolution(SimulationData *pdata);
    ~DiskEvolution();

    void SetInitialCondtion();

    void CalculateDiskGasParamters(double mdot_inf, double r_init, double c_wind);
    void CalculateVelocity();
    void CalculateFlux();

    double CalculateDtMin();
    void CalculateTimeStep(double dt);

private:

    void InterpolationVelocityAtCellBoundary(Array1D<double> &vr_gas_bnd, Array1D<double> &vr_dust_bnd);

    void SetBoundaryCondition();

    void DustGrowth();

    SimulationData *pdata_;
    Reconstruction *prc_;

    Array1D<double> Text_;
    Array1D<double> torque_;
    Array1D<double> Re_;

    Array1D<double> flux_gas_;
    Array1D<double> flux_dust_;
    Array1D<double> flux_dust_m_;

    Array1D<double> mdot_inf_r_;
    Array1D<double> mdot_wind_r_;
    Array1D<double> mdot_r_;

    Array1D<double> dSdt_gas_;
    Array1D<double> dSdt_gas_src_;
    Array1D<double> dSdt_dust_;
    Array1D<double> dSdt_dust_src_;
    Array1D<double> dmSdt_dust_;
    Array1D<double> dmSdt_dust_growth_;

    double mdot_acc_disk_;    // mass accretion rate from disk to star
    double mdot_acc_env_;     // mass accretion rate from envelope to star
    double mdot_acc_total_;   // mdot_acc_disk + mdot_acc_env_
    double mdot_inf_;
    double mdot_wind_;
    double mdot_wind_sweep_;

    double Ls_;
    double Lacc_;


};

#endif /* DISK_EVOLUTION_HPP_ */