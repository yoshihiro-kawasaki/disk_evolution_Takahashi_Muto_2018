#ifndef DRIVER_BASE_HPP_
#define DRIVER_BASE_HPP_

#include "array.hpp"
#include "define_macro.hpp"
#include "sim_data.hpp"
#include "constants.hpp"

class DriverBase
{
public:

    DriverBase(SimulationData *pdata);
    ~DriverBase();

    virtual void SetInitialCondition() = 0;
    virtual void CalculateDiskQuantities() = 0;
    virtual double CalculateDtMin() = 0;
    bool DoTimeStep(const double dt);

protected:

    void CalculateEnclosedMassAndAngularVelocity();
    void CalculateDiskGasQuantities();
    void CalculateGasDustDragCoefficient();
    void CalculateInfallAndWindRate();
    void CalculateGasAndDustVelocity();
    double CalculateDtAdv();
    virtual void CalculateFlux() = 0;
    virtual bool CalculateTimeStep(const double dt) = 0;

    SimulationData *pdata_;

    array::Double1D Text_;
    array::Double1D dOmega_dr_;
    array::Double1D Nvis_;
    array::Double1D dNvis_dr_;
    array::Double1D dP_dr_;
    array::Double1D dPinteg_dr_;
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

    double sigmag_tot_;
    double sigmad_tot_;

private:

};

#endif /* DRIVER_BASE_HPP_ */