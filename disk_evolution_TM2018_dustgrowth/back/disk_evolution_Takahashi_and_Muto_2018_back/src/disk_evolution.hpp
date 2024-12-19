#ifndef DISK_EVOLUTION_HPP_
#define DISK_EVOLUTION_HPP_

#include "array.hpp"
#include "define_macro.hpp"
#include "sim_data.hpp"

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

    SimulationData *pdata_;

    Array1D<double> Text_;
    
    Array1D<double> torque_;

    Array1D<double> flux_gas_;
    Array1D<double> flux_dust_;

    Array1D<double> mdot_inf_r_;
    Array1D<double> mdot_wind_r_;

    Array1D<double> dSdt_gas_;
    Array1D<double> dSdt_dust_;

    double mdot_acc_disk_;    // mass accretion rate from disk to star
    double mdot_acc_env_;     // mass accretion rate from envelope to star
    double mdot_acc_total_;   // mdot_acc_disk + mdot_acc_env_
    double mdot_inf_;
    double mdot_wind_;
    double mdot_wind_sweep_;


};

#endif /* DISK_EVOLUTION_HPP_ */