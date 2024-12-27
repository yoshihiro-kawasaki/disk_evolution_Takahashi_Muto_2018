#include "grid.hpp"

Grid::Grid()
    : is_allocate_arrays_(false)
{

}


Grid::~Grid()
{
    if (is_allocate_arrays_) {
        array::Delete1dArray<double>(r_cen_);
        array::Delete1dArray<double>(r_bnd_);
        array::Delete1dArray<double>(r_vol_);
        array::Delete1dArray<double>(dr_cen_);
        array::Delete1dArray<double>(dr_bnd_);
        array::Delete1dArray<double>(inv_r_cen_);
        array::Delete1dArray<double>(inv_r_vol_);
        array::Delete1dArray<double>(inv_dr_cen_);
        array::Delete2dArray<double>(cdiff_);
        is_allocate_arrays_ = false;
    }
}


void Grid::SetGrid(Utils::InputConfigure &input)
{
    nr_      = input.GetInt("nr");
    rmin_    = input.GetDouble("rmin");
    rmax_    = input.GetDouble("rmax");
    spacing_ = input.GetString("spacing");

    if (is_allocate_arrays_) return;

    ngh_ = 1;
    nrtotc_ = nr_ + 2*ngh_;
    nrtotb_ = nr_ + 2*ngh_ + 1;
    // 計算(disk)領域は、is_ <= n < ie_
    is_     = ngh_;
    ie_     = nr_ + ngh_;

    r_cen_      = array::Allocate1dArray<double>(nrtotc_);
    r_bnd_      = array::Allocate1dArray<double>(nrtotb_);
    r_vol_      = array::Allocate1dArray<double>(nrtotc_);
    dr_cen_     = array::Allocate1dArray<double>(nrtotc_);
    dr_bnd_     = array::Allocate1dArray<double>(nrtotc_);
    inv_r_cen_  = array::Allocate1dArray<double>(nrtotc_);
    inv_r_vol_  = array::Allocate1dArray<double>(nrtotc_);
    inv_dr_cen_ = array::Allocate1dArray<double>(nrtotc_);
    cdiff_      = array::Allocate2dArray<double>(nrtotb_, 3);

    is_allocate_arrays_ = true;

    if (spacing_ == "log") {
        SetLogSpacingGrid();
    } else if (spacing_ == "root") {
        SetRootSpacingGrid();
    }
    SetNumericalDifferentiationCoefficient();
}


void Grid::SetLogSpacingGrid()
{
    double log_rmin = std::log(rmin_);
    double log_rmax = std::log(rmax_);
    double dlogr    = (log_rmax - log_rmin) / (double(nr_));

    for (int i = 0; i < nrtotb_; ++i) {
        r_bnd_[i] = std::exp(log_rmin + dlogr*(double(i) - 1.0)) * cst::ASTRONOMICAL_UNIT;
    }

    for (int i = 0; i < nrtotc_; ++i) {
        r_cen_[i]     = std::sqrt(r_bnd_[i] * r_bnd_[i+1]);
        r_vol_[i]     = M_PI * (SQR(r_bnd_[i+1]) - SQR(r_bnd_[i]));
        inv_r_cen_[i] = 1.0 / r_cen_[i];
        inv_r_vol_[i] = 1.0 / r_vol_[i];
        dr_bnd_[i]    = r_bnd_[i+1] -  r_bnd_[i];
    }

    for (int i = 0; i < nrtotc_-1; ++i) {
        dr_cen_[i]     = r_cen_[i+1] - r_cen_[i];
        inv_dr_cen_[i] = 1.0 / dr_cen_[i];
    }
    return;
}


void Grid::SetRootSpacingGrid()
{
    double p = 0.5;
    double pinv = 1.0 / p;
    double r0p = std::pow(rmin_, p);
    double r1p = std::pow(rmax_, p);
    double drp = (r1p - r0p) / (double(nr_));

    for (int i = 0; i < nrtotb_; ++i) {
        r_bnd_[i] = r0p + drp * i;
        r_bnd_[i] = std::pow(r_bnd_[i], pinv) * cst::ASTRONOMICAL_UNIT;
    }

    for (int i = 0; i < nrtotc_; ++i) {
        r_cen_[i]     = std::sqrt(r_bnd_[i]*r_bnd_[i+1]);
        r_vol_[i]     = M_PI * (SQR(r_bnd_[i+1]) - SQR(r_bnd_[i]));
        inv_r_cen_[i] = 1.0 / r_cen_[i];
        inv_r_vol_[i] = 1.0 / r_vol_[i];
        dr_bnd_[i]    = r_bnd_[i+1] - r_bnd_[i];
    }

    for (int i = 0; i < nrtotc_-1; ++i) {
        dr_cen_[i]     = r_cen_[i+1] - r_cen_[i];
        inv_dr_cen_[i] = 1.0 / dr_cen_[i];
    }

    return;
}

void Grid::CalculateNumericalDifferentiation(const array::Double1D f, array::Double1D dfdr)
{
    for (int i = is_, ie = ie_; i < ie; ++i) {
        dfdr[i] = cdiff_[i][0]*f[i+1] + cdiff_[i][1]*f[i] + cdiff_[i][2]*f[i-1];
    }
    return;
}

void Grid::SetNumericalDifferentiationCoefficient()
{
    int is = is_, ie = ie_;
    double s, t;
    for (int i = is; i < ie; ++i) {
        s = r_cen_[i] - r_cen_[i-1];
        t = r_cen_[i+1] - r_cen_[i];
        cdiff_[i][0] = s*s / (s*t*(s+t));
        cdiff_[i][1] = (t*t - s*s) / (s*t*(s+t));
        cdiff_[i][2] = -t*t / (s*t*(s+t)); 
    }

    return;
}