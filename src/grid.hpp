#ifndef GRID_HPP_
#define GRID_HPP_

#include "array.hpp"
#include "utils.hpp"
#include "constants.hpp"
#include "define_macro.hpp"

#include <iostream>
#include <string>

class Grid 
{
public:

    Grid();
    ~Grid();

    int nr_;                     // number of radial cells
    int ngh_;                    // number of ghost cells
    int nrtotc_;                 // total cell center
    int nrtotb_;                 // total cell boundary
    int is_;                     // start index of computation domain
    int ie_;                     // end index of computation domain
    double rmin_;                // minimum radius
    double rmax_;                // maximum radius 
    array::Double1D r_cen_;      // cell center [cm]
    array::Double1D r_bnd_;      // cell boundary [cm]
    array::Double1D r_vol_;      // cell volume [cm]
    array::Double1D dr_cen_;     // delta cell center [cm]
    array::Double1D dr_bnd_;     // delta cell boundary [cm]
    array::Double1D inv_r_cen_;  // 1 / r_cen_
    array::Double1D inv_r_vol_;  // 1 / r_vol_
    array::Double1D inv_dr_cen_; // 1 / dr_cen_

    void SetGrid(Utils::InputConfigure &input);
    void CalculateNumericalDifferentiation(const array::Double1D f, array::Double1D dfdr);

private:

    array::Double2D cdiff_;      // coefficient for numerical derivative
    std::string spacing_;
    bool is_allocate_arrays_;

    void SetLogSpacingGrid();
    void SetRootSpacingGrid();
    void SetNumericalDifferentiationCoefficient();
};

#endif /* GRID_HPP_ */