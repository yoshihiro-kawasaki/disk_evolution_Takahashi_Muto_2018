#ifndef DRIVER_HPP_
#define DRIVER_HPP_


#include "driver_base.hpp"

class Driver
    : public DriverBase
{
public:
    Driver(SimulationData *pdata);
    ~Driver();
    
    void SetInitialCondition();
    double CalculateDtMin();
    void CalculateDiskQuantities();

protected:

    void CalculateDustStokesNumber();
    void CalculateDustGrowthRate();
    void CalculateFlux();
    bool CalculateTimeStep(const double dt);

};

#endif /* DRIVER_HPP_ */