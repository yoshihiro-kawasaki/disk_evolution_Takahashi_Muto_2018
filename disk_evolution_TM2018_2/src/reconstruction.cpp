#include "reconstruction.hpp"

Reconstruction::Reconstruction(SimulationData *pdata)
    : pdata_(pdata)
{

}

Reconstruction::~Reconstruction()
{

}

void Reconstruction::VanLeerReconstruction(const Array1D<double> &Q, Array1D<double> &Ql, Array1D<double> &Qr)
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;

    double Qpp, Qpm, dQdr_mono;

    for (int i = is; i <= ie; ++i) {

        Qpm = (Q(i-1) - Q(i-2)) * pdata_->grid_.inv_dr_cen_(i-2);
        Qpp = (Q(i) - Q(i-1)) * pdata_->grid_.inv_dr_cen_(i-1);
        if (Qpm + Qpp == 0.0) {
            dQdr_mono = 0.0;
        } else {
            dQdr_mono = 2.0 * Qpm * Qpp / (Qpm + Qpp);
        }
        Ql(i) = Q(i-1) + (pdata_->grid_.r_bnd_(i) - pdata_->grid_.r_cen_(i-1)) * dQdr_mono;

        Qpm = (Q(i) - Q(i-1)) * pdata_->grid_.inv_dr_cen_(i-1);
        Qpp = (Q(i+1) - Q(i)) * pdata_->grid_.inv_dr_cen_(i);
        dQdr_mono = 2.0 * Qpm * Qpp / (Qpm + Qpp);
        if (Qpm + Qpp == 0.0) {
            dQdr_mono = 0.0;
        } else {
            dQdr_mono = 2.0 * Qpm * Qpp / (Qpm + Qpp);
        }
        Qr(i) = Q(i) + (pdata_->grid_.r_bnd_(i) - pdata_->grid_.r_cen_(i)) * dQdr_mono;
    }

    return;
}