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
        pdata_->dust_.md_[i]          = 0.0;
        pdata_->dust_.ad_[i]          = 0.0;
        pdata_->dust_.mSigma_[i]      = 0.0;
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
    CalculateDustGrowthRate();
    return;
}


void Driver::CalculateDustStokesNumber()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;
    double inv_mfp;

    for (int i = is; i < ie; ++i) {
        if (pdata_->gas_.Sigma_[i] < SIGMA_GAS_FLOOR) {
            pdata_->dust_.St_[i] = 0.0;
            continue;
        }
        inv_mfp                  = (4.0 * pdata_->gas_.rho_mid_[i] * cst::H2_CROSS_SECTION) / (9.0 * cst::GAS_MOLECULAR_MASS);
        pdata_->dust_.St_[i] = (0.5 * M_PI * pdata_->dust_.rho_di_ * pdata_->dust_.ad_[i] / pdata_->gas_.Sigma_[i]) * std::max(1.0, pdata_->dust_.ad_[i]*inv_mfp);
    }
    return;
}


void Driver::CalculateDustGrowthRate()
{
    int is = pdata_->grid_.is_,
        ie = pdata_->grid_.ie_;

    double StP, StT, HdT, Re, Stmin, StTPast, sqr_dvbr, sqr_dvturbI, sqr_dvturbII, sqr_dvturb, dv, deltam;

    for (int i = is; i < ie; ++i) {
        if (pdata_->gas_.Sigma_[i] < SIGMA_GAS_FLOOR || pdata_->dust_.Sigma_[i] < SIGMA_DUST_FLOOR) {
            pdata_->dust_.Hd_[i]  = 0.0;
            dmSdt_dust_growth_[i] = 0.0;
            continue;
        }

        StP = 0.5 * pdata_->dust_.St_[i];
        StT = pdata_->dust_.St_[i];

        pdata_->dust_.Hd_[i] = pdata_->gas_.Hg_[i] / std::sqrt(1.0 + StT*(1.0 + 2.0*StT) / (pdata_->gas_.alpha_turb_ * (1.0 + StT)));
        HdT                  = pdata_->dust_.Hd_[i];

        sqr_dvbr     = 2.0 * cst::COE_VBR * pdata_->gas_.T_[i] / pdata_->dust_.md_[i];

        Re           = 0.5 * pdata_->gas_.alpha_turb_ * pdata_->gas_.Sigma_[i] * cst::H2_CROSS_SECTION * cst::INV_MG;
        Stmin        = 1.0 / Re;
        StTPast      = std::max(Stmin, std::min(1.6*StT, 1.0));
        sqr_dvturbI  = (StT-StP)/(StT + StP) * (SQR(StT)/(StTPast + StT) - SQR(StT)/(1+StT) - SQR(StP)/(StTPast + StP) + SQR(StP)/(1+StP));
        sqr_dvturbII = (StTPast - Stmin) * (((StTPast + Stmin)*StT + StTPast*Stmin)/((StT+StTPast)*(StT+Stmin)) + ((StTPast + Stmin)*StP + StTPast*Stmin)/((StP+StTPast)*(StP+Stmin)));
        sqr_dvturb   = pdata_->gas_.alpha_turb_ * SQR(pdata_->gas_.cs_[i]) * (sqr_dvturbI + sqr_dvturbII);
        if (sqr_dvturb < 0.0) sqr_dvturb = 0.0;

        dv = std::sqrt(sqr_dvbr + sqr_dvturb);

        deltam = std::min(1.0, - std::log(dv/pdata_->dust_.vfrag_)/std::log(5));

        dmSdt_dust_growth_[i] = 2.0 * cst::SQRT_PI * SQR(pdata_->dust_.ad_[i]) * dv * SQR(pdata_->dust_.Sigma_[i]) * deltam / HdT;
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
                flux_dust_[i]   = 0.0;
                flux_dust_m_[i] = 0.0;
            } else {
                flux_dust_[i]   = pdata_->dust_.Sigma_[i-1]  * vdr_bnd;
                flux_dust_m_[i] = pdata_->dust_.mSigma_[i-1] * vdr_bnd;  
            }
        } else {
            if (i == ie) {
                flux_dust_[i]   = 0.0;
                flux_dust_m_[i] = 0.0;
            } else {
                flux_dust_[i]   = pdata_->dust_.Sigma_[i]  * vdr_bnd;
                flux_dust_m_[i] = pdata_->dust_.mSigma_[i] * vdr_bnd;
            }
        }

    }

    return;
}

double Driver::CalculateDtMin() 
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;
    double dt_temp, inv_dV, dSdt_gas_adv, dSdt_gas_src, dSdt_dust_adv, dSdt_dust_src, dmSdt_dust_adv, dmSdt_dust_src, dmSdt_dust_growth, m;

    CalculateFlux();

    dt_temp  = 1.0e100;

    for (int i = is; i < ie; ++i) {

        inv_dV = pdata_->grid_.inv_r_vol_[i];

        // gas
        dSdt_gas_adv = - 2.0*M_PI * (pdata_->grid_.r_bnd_[i+1]*flux_gas_[i+1] - pdata_->grid_.r_bnd_[i]*flux_gas_[i]) * inv_dV;
        dSdt_gas_src = ((1.0 - pdata_->dust_.fdg_in_) * (mdot_inf_r_[i+1] - mdot_inf_r_[i]) - (mdot_wind_r_[i+1] - mdot_wind_r_[i])) * inv_dV;
        dSdt_gas_[i] = dSdt_gas_adv + dSdt_gas_src;

        if (pdata_->gas_.vr_[i] != 0.0 && pdata_->gas_.Sigma_[i] > SIGMA_GAS_FLOOR) {
            dt_temp = std::min(dt_temp, pdata_->grid_.dr_bnd_[i] / std::abs(pdata_->gas_.vr_[i]));
        }
        dt_temp = std::min(dt_temp, 1.0/pdata_->gas_.Omega_[i]);

        if (dSdt_gas_src != 0.0 && pdata_->gas_.Sigma_[i] > SIGMA_GAS_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->gas_.Sigma_[i]/dSdt_gas_src));
        }

        // dust
        dSdt_dust_adv = - 2.0*M_PI* (pdata_->grid_.r_bnd_[i+1]*flux_dust_[i+1] - pdata_->grid_.r_bnd_[i]*flux_dust_[i]) * inv_dV;
        dSdt_dust_src = pdata_->dust_.fdg_in_*(mdot_inf_r_[i+1] - mdot_inf_r_[i]) * inv_dV;
        dSdt_dust_[i] = dSdt_dust_adv + dSdt_dust_src;

        dmSdt_dust_adv    = - 2.0*M_PI* (pdata_->grid_.r_bnd_[i+1]*flux_dust_m_[i+1] - pdata_->grid_.r_bnd_[i]*flux_dust_m_[i]) * inv_dV;
        dmSdt_dust_growth = dmSdt_dust_growth_[i];
        if (pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            m = pdata_->dust_.md_[i];
        } else {
            m = pdata_->dust_.md0_;
        }
        dmSdt_dust_src = m * dSdt_dust_src;
        dmSdt_dust_[i] = dmSdt_dust_adv + dmSdt_dust_growth + dmSdt_dust_src;

        if (pdata_->dust_.vr_[i] != 0.0 && pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->grid_.dr_bnd_[i] / pdata_->dust_.vr_[i]));
        }

        if (dSdt_dust_src != 0.0 && pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->dust_.Sigma_[i]/dSdt_dust_src));
        }

        if (dmSdt_dust_growth != 0.0 & pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->dust_.mSigma_[i] / dmSdt_dust_growth));
        }
    }

    mdot_acc_disk_  = - 2.0*M_PI*pdata_->grid_.r_bnd_[is]*(flux_gas_[is] + flux_dust_[is]);
    mdot_acc_env_   = mdot_inf_r_[is];
    mdot_acc_total_ = mdot_acc_disk_ + mdot_acc_env_;

    return dt_temp;
}


bool Driver::CalculateTimeStep(const double dt)
{
    int is = pdata_->grid_.is_, 
        ie = pdata_->grid_.ie_;
    double dV, inv_rho_di, one_third;


    sigmag_tot_ = 0.0;
    sigmad_tot_ = 0.0;

    inv_rho_di = 1.0 / pdata_->dust_.rho_di_;
    one_third  = 1.0 / 3.0;

    for (int i = is; i < ie; ++i) {
        pdata_->gas_.Sigma_[i]   += dSdt_gas_[i] * dt;
        pdata_->dust_.Sigma_[i]  += dSdt_dust_[i] * dt;
        pdata_->dust_.mSigma_[i] += dmSdt_dust_[i] * dt;
        if (pdata_->dust_.Sigma_[i] > SIGMA_DUST_FLOOR) {
            pdata_->dust_.md_[i] = pdata_->dust_.mSigma_[i] / pdata_->dust_.Sigma_[i];
            pdata_->dust_.ad_[i] = std::pow(cst::COE_INV_SPHERE * pdata_->dust_.md_[i] * inv_rho_di, one_third);
        }

        dV           = pdata_->grid_.r_vol_[i];
        sigmag_tot_ += pdata_->gas_.Sigma_[i] * dV;
        sigmad_tot_ += pdata_->dust_.Sigma_[i] * dV;
    }

    return true;
}