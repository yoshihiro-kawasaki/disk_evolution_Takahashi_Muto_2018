#include "../driver.hpp"

Driver::Driver(SimulationData *pdata)
    : DriverBase(pdata)
{
    
}


Driver::~Driver()
{

}


void Driver::SetInitialCondition()
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;
    double temp;
    for (int i = is; i < ie; ++i) {
        pdata_->gas_.Sigma_[i]        = 0.0;
        pdata_->gas_.Sigma_floor_[i]  = SIGMA_GAS_FLOOR;
        pdata_->dust_.Sigma_[i]       = 0.0;
        pdata_->dust_.Sigma_floor_[i] = SIGMA_DUST_FLOOR;
        temp                          = 1.5e2 * std::pow(pdata_->grid_.r_cen_[i]/cst::ASTRONOMICAL_UNIT, -0.42857142857142855);
        Text_[i]                      = std::max(temp, 10.0);
    }
    return;
}


void Driver::CalculateDiskQuantities()
{
    CalculateEnclosedMassAndAngularVelocity();
    CalculateDiskGasQuantities();
    CalculateDustStokesNumber();
    CalculateGasDustDragCoefficient();
    CalculateInfallAndWindRate();
    CalculateGasAndDustVelocity();
    return;
}


void Driver::CalculateDustStokesNumber()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;
    double mfp, cSt;

    cSt = 0.5 * M_PI * pdata_->dust_.rho_di_ * pdata_->dust_.ad0_;

    for (int i = is; i < ie; ++i) {
        if (pdata_->gas_.Sigma_[i] < SIGMA_GAS_FLOOR) {
            pdata_->dust_.St_[i] = 0.0;
            continue;
        }

        mfp = cst::GAS_MOLECULAR_MASS / (pdata_->gas_.rho_mid_[i] * cst::H2_CROSS_SECTION);
        pdata_->dust_.St_[i] = (cSt / pdata_->gas_.Sigma_[i]) * std::max(1.0, 4.0*pdata_->dust_.ad0_/(9.0*mfp));
    }
    
    return;
}


void Driver::CalculateFlux()
{
    int is     = pdata_->grid_.is_, 
        ie     = pdata_->grid_.ie_,
        nrtotb = pdata_->grid_.nrtotb_;
    double rb_rc, m, vgr_bnd, vdr_bnd;

    // zero gradient at boudaries
    pdata_->gas_.vr_[is-1]  = pdata_->gas_.vr_[is]; 
    pdata_->dust_.vr_[is-1] = pdata_->dust_.vr_[is];
    pdata_->gas_.vr_[ie]    = pdata_->gas_.vr_[ie-1];
    pdata_->dust_.vr_[ie]   = pdata_->dust_.vr_[ie-1];

    for (int i = is; i <= ie; ++i) {

        rb_rc = pdata_->grid_.r_bnd_[i] - pdata_->grid_.r_cen_[i-1];

        // gas
        m       = (pdata_->gas_.vr_[i] -  pdata_->gas_.vr_[i-1]) * pdata_->grid_.inv_dr_cen_[i-1];
        vgr_bnd = pdata_->gas_.vr_[i-1] + m*rb_rc;

        if (vgr_bnd >= 0.0) {
            if (i == is) {
                flux_gas_[i] = 0.0;
            } else {
                flux_gas_[i] = pdata_->gas_.Sigma_[i-1] * vgr_bnd;
            }
        } else {
            if (i == ie) {
                flux_gas_[i] = 0.0;
            } else {
                flux_gas_[i] = pdata_->gas_.Sigma_[i] * vgr_bnd;
            }
        }

        // dust
        m       = (pdata_->dust_.vr_[i] -  pdata_->dust_.vr_[i-1]) * pdata_->grid_.inv_dr_cen_[i-1];
        vdr_bnd = pdata_->dust_.vr_[i-1] + m*rb_rc;

        if (vdr_bnd >= 0.0) {
            if (i == is) {
                flux_dust_[i] = 0.0;
            } else {
                flux_dust_[i] = pdata_->dust_.Sigma_[i-1] * vdr_bnd;
            }
        } else {
            if (i == ie) {
                flux_dust_[i] = 0.0;
            } else {
                flux_dust_[i] = pdata_->dust_.Sigma_[i] * vdr_bnd;
            }
        }

    }

    return;
}

double Driver::CalculateDtMin() 
{ 
    return CalculateDtAdv(); 
}


bool Driver::CalculateTimeStep(const double dt)
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;
    double dV;

    sigmag_tot_ = 0.0;
    sigmad_tot_ = 0.0;

    for (int i = is; i < ie; ++i) {
        pdata_->gas_.Sigma_[i]  += dSdt_gas_[i] * dt;
        pdata_->dust_.Sigma_[i] += dSdt_dust_[i] * dt;

        dV           = pdata_->grid_.r_vol_[i];
        sigmag_tot_ += pdata_->gas_.Sigma_[i] * dV;
        sigmad_tot_ += pdata_->dust_.Sigma_[i] * dV;
    }

    return true;
}