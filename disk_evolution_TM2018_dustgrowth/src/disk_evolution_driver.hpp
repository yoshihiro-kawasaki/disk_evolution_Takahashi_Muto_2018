#ifndef DISK_EVOLUTION_HPP_
#define DISK_EVOLUTION_HPP_

#include "array.hpp"
#include "define_macro.hpp"
#include "sim_data.hpp"
#include "constants.hpp"

constexpr double COE_SPHERE     = 4.0 * M_PI / 3.0;
constexpr double COE_INV_SPHERE = 1.0 / COE_SPHERE;
constexpr double COE_VBR        = 8.0 * cst::BOLTZMANN_CONSTANT / M_PI;

class DiskEvolutionDriver
{
public:

    DiskEvolutionDriver(SimulationData *pdata);
    ~DiskEvolutionDriver();

    void SetInitialCondtion();
    void CalculateDiskQuantities();
    double CalculateDtMin();
    bool CalculateTimeStep(const double dt);

private:

    void CalculateEnclosedMassAndAngularVelocity();
    void CalculateDiskGasQuantities();
    void CalculateGasDustDragCoefficient();
    void CalculateDustGrowthRate();
    void CalculateInfallAndWindRate();
    void CalculateGasAndDustVelocity();
    void CalculateFlux();

    SimulationData *pdata_;

    array::Double1D Text_;
    array::Double1D dOmega_dr_;
    array::Double1D Nvis_;
    array::Double1D dNvis_dr_;
    array::Double1D dP_dr_;
    array::Double1D dPinteg_dr_;
    array::Double1D Re_;
    array::Double1D flux_gas_;
    array::Double1D flux_dust_;
    array::Double1D flux_dust_m_;
    array::Double1D mdot_inf_r_;
    array::Double1D mdot_wind_r_;
    array::Double1D mdot_r_;
    array::Double1D dSdt_gas_;
    array::Double1D dSdt_gas_src_;
    array::Double1D dSdt_dust_;
    array::Double1D dSdt_dust_src_;
    array::Double1D dmSdt_dust_;
    array::Double1D dmSdt_dust_growth_;

    bool is_allocate_arrays_;

    double mdot_acc_disk_;    // mass accretion rate from disk to star
    double mdot_acc_env_;     // mass accretion rate from envelope to star
    double mdot_acc_total_;   // mdot_acc_disk + mdot_acc_env_
    double mdot_inf_;         // infall rate from envelope to disk
    double mdot_wind_;        // wind mass loss rate from disk
};

#endif /* DISK_EVOLUTION_HPP_ */