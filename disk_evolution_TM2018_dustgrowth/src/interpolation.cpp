#include "interpolation.hpp"

Interpolation::Interpolation(SimulationData *pdata)
    : pdata_(pdata)
{
    cp_.Resize(4);

    int is = pdata_->grid_.is_;

    cp_(1) = (pdata_->grid_.r_cen_(is-1) - pdata_->grid_.r_cen_(is+1)) * (pdata_->grid_.r_cen_(is-1) - pdata_->grid_.r_cen_(is+2))
        / ((pdata_->grid_.r_cen_(is) - pdata_->grid_.r_cen_(is+1)) * (pdata_->grid_.r_cen_(is) - pdata_->grid_.r_cen_(is+2)));

    cp_(2) = (pdata_->grid_.r_cen_(is-1) - pdata_->grid_.r_cen_(is)) * (pdata_->grid_.r_cen_(is-1) - pdata_->grid_.r_cen_(is+2))
        / ((pdata_->grid_.r_cen_(is+1) - pdata_->grid_.r_cen_(is)) * (pdata_->grid_.r_cen_(is+1) - pdata_->grid_.r_cen_(is+2)));

    cp_(3) = (pdata_->grid_.r_cen_(is-1) - pdata_->grid_.r_cen_(is)) * (pdata_->grid_.r_cen_(is-1) - pdata_->grid_.r_cen_(is+1))
        / ((pdata_->grid_.r_cen_(is+2) - pdata_->grid_.r_cen_(is)) * (pdata_->grid_.r_cen_(is+2) - pdata_->grid_.r_cen_(is+1)));

}