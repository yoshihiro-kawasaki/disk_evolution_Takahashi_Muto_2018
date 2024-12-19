#ifndef RECONSTRUCTION_HPP_
#define RECONSTRUCTION_HPP_

#include "array.hpp"
#include "sim_data.hpp"

class Reconstruction
{
public:
    Reconstruction(SimulationData *pdata);
    ~Reconstruction();

    void VanLeerReconstruction(const Array1D<double> &Q, Array1D<double> &Ql, Array1D<double> &Qr);

private:
    SimulationData *pdata_;
};

#endif /* RECONSTRUCTION_HPP_ */