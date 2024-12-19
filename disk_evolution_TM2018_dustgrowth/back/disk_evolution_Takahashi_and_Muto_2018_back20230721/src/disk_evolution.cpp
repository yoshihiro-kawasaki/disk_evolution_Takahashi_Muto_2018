#include "disk_evolution.hpp"

DiskEvolution::DiskEvolution(SimulationData *pdata)
    : pdata_(pdata)
{
    int nrtotc = pdata_->grid_.nrtotc_;
    int nrtotb = pdata_->grid_.nrtotb_;

    Text_.Resize(nrtotc);

    torque_.Resize(nrtotc);

    cvg_.Resize(nrtotc, 3);
    cvd_.Resize(nrtotb, 4);

    flux_gas_.Resize(nrtotb);
    flux_dust_.Resize(nrtotb);

    mdot_inf_r_.Resize(nrtotb);
    mdot_wind_r_.Resize(nrtotb);

    dSdt_gas_.Resize(nrtotc);
    dSdt_dust_.Resize(nrtotc);

    mdotloss_ = 0.0;
}

DiskEvolution::~DiskEvolution()
{

}


void DiskEvolution::SetInitialCondtion()
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;

    double temp;

    for (int i = is; i < ie; ++i) {
        pdata_->gas_.Sigma_(i)       = 0.0;
        pdata_->gas_.Sigma_floor_(i) = 1.0e-20;
        pdata_->dust_.Sigma_(i)      = 0.0;
        temp                         = 1.5e2 * std::pow(pdata_->grid_.r_cen_(i)*inv_AU, -0.42857142857142855);
        Text_(i) = std::max(temp, 10.0);
    }

    return;
}


void DiskEvolution::CalculateDiskGasParamters(double mdot_inf, double r_init, double c_wind)
{
    /*
    円盤内の物理量を計算

    Parameters
    mdot_inf : total infalling gas mass rate
    r_init   : initial cloud shell radius 
    c_wind   : efficiency of disk wind

    */

    int is = pdata_->grid_.is_;
    int ie = pdata_->grid_.ie_;

    double encl_mr, inv_encl_mr, dOmega_dr, rhog_mid, temp, kappa_epy, Qt, inv_Sigma_gas, St, nu;
    double fmass, coe1, coe2, lambda1, lambda2;

    // enclosed mass at cell boundary
    pdata_->gas_.Mr_(is) = pdata_->star_.mass_;
    for (int i = is+1; i <= ie; ++i) {
        pdata_->gas_.Mr_(i) = pdata_->gas_.Mr_(i-1) + pdata_->grid_.r_vol_(i-1) * pdata_->gas_.Sigma_(i-1);
    }

    // others
    for (int i = is; i < ie; ++i) {

        // enclosed mass at cell center 
        encl_mr = std::sqrt(pdata_->gas_.Mr_(i+1) * pdata_->gas_.Mr_(i));
        // encl_mr = 0.5 * (pdata_->gas_.Mr_(i+1) + pdata_->gas_.Mr_(i));
        inv_encl_mr = 1.0 / encl_mr;

        pdata_->gas_.Omega_(i) = std::sqrt(pc::GRAVITATIONAL_CONSTANT*encl_mr * CUB(pdata_->grid_.inv_r_cen_(i)));
        dOmega_dr = pdata_->gas_.Omega_(i) * (M_PI*pdata_->grid_.r_cen_(i)*pdata_->gas_.Sigma_(i) * inv_encl_mr - 1.5*pdata_->grid_.inv_r_cen_(i));
        pdata_->gas_.j_(i) = SQR(pdata_->grid_.r_cen_(i)) * pdata_->gas_.Omega_(i);
        
        if (pdata_->gas_.Sigma_(i) < pdata_->gas_.Sigma_floor_(i)) {

            pdata_->gas_.T_(i)     = 10.0;
            pdata_->gas_.cs_(i)    = std::sqrt(pc::BOLTZMANN_CONSTANT * 10.0 * inv_mg);
            pdata_->gas_.Hg_(i)    = pdata_->gas_.cs_(i) / pdata_->gas_.Omega_(i);
            rhog_mid               = 0.0;
            pdata_->gas_.P_(i)     = 0.0;
            Qt                     = 1e3;
            pdata_->gas_.alpha_(i) = 0.0;
            torque_(i)             = 0.0;

            cvg_(i, 0) = 0.0;
            cvg_(i, 1) = 0.0;
            cvg_(i, 2) = 0.0;

            cvd_(i, 0) = 0.0;
            cvd_(i, 1) = 0.0;
            cvd_(i, 2) = 0.0;
            cvd_(i, 3) = 0.0;


        } else {

            pdata_->gas_.T_(i) = Text_(i);

            pdata_->gas_.cs_(i) = std::sqrt(pc::BOLTZMANN_CONSTANT * pdata_->gas_.T_(i) * inv_mg);

            pdata_->gas_.Hg_(i) = pdata_->gas_.cs_(i) / pdata_->gas_.Omega_(i);

            rhog_mid = pdata_->gas_.Sigma_(i) / (sqrt_2pi * pdata_->gas_.Hg_(i));

            pdata_->gas_.P_(i) = SQR(pdata_->gas_.cs_(i)) * rhog_mid;

            // kappa_epy = pdata_->gas_.Omega_(i)*std::sqrt(1.0 + 2.0*M_PI*SQR(pdata_->grid_.r_cen_(i))*pdata_->gas_.Sigma_(i) * inv_encl_mr);
            // Qt = kappa_epy * pdata_->gas_.cs_(i) / (M_PI * pc::GRAVITATIONAL_CONSTANT * pdata_->gas_.Sigma_(i));
            Qt = pdata_->gas_.Omega_(i) * pdata_->gas_.cs_(i) / (M_PI * pc::GRAVITATIONAL_CONSTANT * pdata_->gas_.Sigma_(i));

            pdata_->gas_.alpha_(i) = std::exp(-SQR(Qt)*SQR(Qt)) + pdata_->gas_.alpha_turb_;

            // torque_(i) = - SQR(pdata_->grid_.r_cen_(i)) * pdata_->gas_.alpha_(i) * SQR(pdata_->gas_.cs_(i)) * pdata_->gas_.Sigma_(i);

            nu = pdata_->gas_.alpha_(i) * pdata_->gas_.cs_(i) * pdata_->gas_.Hg_(i);
            torque_(i) = CUB(pdata_->grid_.r_cen_(i)) * nu * pdata_->gas_.Sigma_(i) * dOmega_dr;

            inv_Sigma_gas = 1.0 / pdata_->gas_.Sigma_(i);
            
            St = pdata_->dust_.St_;

            fmass = 2.0 * M_PI * SQR(pdata_->grid_.r_cen_(i)) * pdata_->gas_.Sigma_(i) * inv_encl_mr;
            coe1 = 1.0 / (1.0 + (1.0 + fmass)*SQR(St));

            cvd_(i, 0) = (1.0 + fmass*SQR(St)) * coe1;
            cvd_(i, 1) = St * coe1;
            cvd_(i, 2) = SQR(St) * coe1;
            cvd_(i, 3) = coe1;

            lambda1 = St * pdata_->dust_.Sigma_(i) * coe1 * inv_Sigma_gas;
            lambda2 = pdata_->dust_.Sigma_(i) * coe1 * inv_Sigma_gas;
            coe2 = 1.0 / (SQR(1.0+lambda2) + SQR(lambda1)*(1.0 + fmass)); 

            cvg_(i, 0) = (1.0 + lambda2) * coe2;
            cvg_(i, 1) = 2.0*lambda1 * (1.0 + fmass) * coe2;
            cvg_(i, 2) = lambda1 * coe2 * 0.5;

        } // end if

    } // end for


    // radial velocity at cell center due to viscous torque and pressure gradient
    double gradT, gradP, inv_dr, eta;

    // 境界条件(ゴーストセルでの値) : 改善の余地あり
    torque_(0) = torque_(1); // zero torque at inner boundary
    // torque_(0) = torque_(1) 
    //     + (torque_(2) - torque_(1)) / (pdata_->grid_.r_cen_(2) - pdata_->grid_.r_cen_(1)) 
    //     * (pdata_->grid_.r_cen_(0) - pdata_->grid_.r_cen_(1));
    // torque_(ie) = torque_(ie-2);
    pdata_->gas_.P_(0)  = pdata_->gas_.P_(1); 
        + (pdata_->gas_.P_(2) - pdata_->gas_.P_(1)) / (pdata_->grid_.r_cen_(2) - pdata_->grid_.r_cen_(1)) 
        * (pdata_->grid_.r_cen_(0) - pdata_->grid_.r_cen_(1));
    pdata_->gas_.P_(ie) = pdata_->gas_.P_(ie-2);

    for (int i = is; i < ie; ++i) {

        if (pdata_->gas_.Sigma_(i) < pdata_->gas_.Sigma_floor_(i)) {
            pdata_->gas_.vr_vis_(i) = 0.0;
            pdata_->gas_.vr_eta_(i) = 0.0;
            continue;
        }

        inv_dr = 1.0 / (pdata_->grid_.r_cen_(i+1) -  pdata_->grid_.r_cen_(i-1));
        gradT = (torque_(i+1) - torque_(i-1)) * inv_dr;
        gradP = (pdata_->grid_.r_cen_(i) /  pdata_->gas_.P_(i)) * (pdata_->gas_.P_(i+1) - pdata_->gas_.P_(i-1)) * inv_dr;

        pdata_->gas_.vr_vis_(i) = 2.0 * gradT / (pdata_->gas_.j_(i) * pdata_->gas_.Sigma_(i));
        eta = - 0.5 * SQR(pdata_->gas_.Hg_(i) * pdata_->grid_.inv_r_cen_(i)) * gradP;
        pdata_->gas_.vr_eta_(i)  = eta * pdata_->grid_.r_cen_(i) * pdata_->gas_.Omega_(i);

    }


    /* calculate mdot_inf and flux */
    double j, j_core, inv_j_core;
    j_core = pdata_->cloud_.Omega_c_ * SQR(r_init);
    inv_j_core = 1.0 / j_core;

    for (int i = is; i <= ie; ++i) {
        j = std::sqrt(pc::GRAVITATIONAL_CONSTANT*pdata_->gas_.Mr_(i)*pdata_->grid_.r_bnd_(i));
        if (j <= j_core) {
            mdot_inf_r_(i) = mdot_inf * (1.0 - std::sqrt(1.0 - j*inv_j_core));
        } else {
            mdot_inf_r_(i) = mdot_inf;
        }
    }

    // disk wind  Takahashi and Muto 2018
    if (c_wind != 0.0) {
        mdot_wind_r_(1) = 0.0;
        for (int i = is+1; i <= ie; ++i) {
            mdot_wind_r_(i) = mdot_wind_r_(i-1) + pdata_->grid_.r_vol_(i-1) * c_wind * pdata_->gas_.Sigma_(i-1) * pdata_->gas_.Omega_(i-1);
        }
    }

    double mdot_tot_r, vgphi;
    for (int i = is; i < ie; ++i) {

        // velocity due to infall or wind
        mdot_tot_r = 0.5 * (mdot_inf_r_(i) + mdot_inf_r_(i+1) - (mdot_wind_r_(i)  + mdot_wind_r_(i+1)));
        encl_mr = 0.5*(pdata_->gas_.Mr_(i) + pdata_->gas_.Mr_(i+1));
        pdata_->gas_.vr_src_(i)  = - pdata_->grid_.r_cen_(i) * mdot_tot_r / encl_mr;

        // total radial velocity
        pdata_->gas_.vr_(i)  = cvg_(i, 0)*pdata_->gas_.vr_vis_(i)  
                             + cvg_(i, 1)*pdata_->gas_.vr_eta_(i)  
                             + pdata_->gas_.vr_src_(i);
        vgphi = cvg_(i, 2)*pdata_->gas_.vr_vis_(i) 
              - cvg_(i, 0)*pdata_->gas_.vr_eta_(i);

        pdata_->dust_.vr_(i) = cvd_(i, 0)*pdata_->gas_.vr_(i) 
                             + 2.0*cvd_(i, 1)*vgphi 
                             + cvd_(i, 2)*pdata_->gas_.vr_src_(i);
    }

    return;

}

void DiskEvolution::InterpolationVelocityAtCellBoundary(Array1D<double> &vr_gas_bnd, Array1D<double> &vr_dust_bnd)
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;
    double rc, rb, inv_dr, m;

    // i = is
    rc = pdata_->grid_.r_cen_(is);
    rb = pdata_->grid_.r_bnd_(is);
    inv_dr = 1.0 / (pdata_->grid_.r_cen_(is+1) - pdata_->grid_.r_cen_(is));
    // gas
    m = (pdata_->gas_.vr_(is+1) - pdata_->gas_.vr_(is)) * inv_dr;
    vr_gas_bnd(is)  = pdata_->gas_.vr_(is) + m*(rb - rc);
    // dust
    vr_dust_bnd(is) = pdata_->dust_.vr_(is) + m*(rb - rc);

    vr_gas_bnd(ie)  = 0.0;
    vr_dust_bnd(ie) = 0.0;


    for (int i = is+1; i < ie; ++i) {
        rc = pdata_->grid_.r_cen_(i-1);
        rb = pdata_->grid_.r_bnd_(i);
        inv_dr = 1.0 / (pdata_->grid_.r_cen_(i) -  pdata_->grid_.r_cen_(i-1));

        // gas 
        m = (pdata_->gas_.vr_(i) -  pdata_->gas_.vr_(i-1)) * inv_dr;
        vr_gas_bnd(i) = pdata_->gas_.vr_(i-1) + m*(rb - rc);

        // dust
        m = (pdata_->dust_.vr_(i) -  pdata_->dust_.vr_(i-1)) * inv_dr;
        vr_dust_bnd(i) = pdata_->dust_.vr_(i-1) + m*(rb - rc);
    }

    return;
}

void DiskEvolution::CalculateFlux()
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;
    Array1D<double> vr_gas_bnd(pdata_->grid_.nrtotb_), vr_dust_bnd(pdata_->grid_.nrtotb_);

    InterpolationVelocityAtCellBoundary(vr_gas_bnd, vr_dust_bnd);

    for (int i = is; i <= ie; ++i) {

        // gas
        if (vr_gas_bnd(i) >= 0.0) {
            if (i == is) {
                flux_gas_(i) = 0.0;
            } else {
                flux_gas_(i) = pdata_->gas_.Sigma_(i-1) * vr_gas_bnd(i);
            }
        } else {
            if (i == ie) {
                flux_gas_(i) == 0.0;
            } else {
                flux_gas_(i) = pdata_->gas_.Sigma_(i) * vr_gas_bnd(i);
            }
        }

        // dust
        if (vr_dust_bnd(i) >= 0.0) {
            if (i == is) {
                flux_dust_(i) = 0.0;
            } else {
                flux_dust_(i) = pdata_->dust_.Sigma_(i-1) * vr_dust_bnd(i);
            }
        } else {
            if (i == ie) {
                flux_dust_(i) == 0.0;
            } else {
                flux_dust_(i) = pdata_->dust_.Sigma_(i) * vr_dust_bnd(i);
            }
        }
    }

    return;
}


double DiskEvolution::CalculateDtMin()
{
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;
    double dt_temp, dV_inv, dSdt_gas_adv, dSdt_gas_src, dSdt_dust_adv, dSdt_dust_src;

    dt_temp  = 1.0e100;

    for (int i = is; i < ie; ++i) {

        dV_inv = pdata_->grid_.inv_r_vol_(i);

        // gas

        dSdt_gas_adv  = - 2.0*M_PI * (pdata_->grid_.r_bnd_(i+1)*flux_gas_(i+1) - pdata_->grid_.r_bnd_(i)*flux_gas_(i)) * dV_inv;
        dSdt_gas_src = (mdot_inf_r_(i+1) - mdot_inf_r_(i) - (mdot_wind_r_(i+1) - mdot_wind_r_(i))) * dV_inv;
        dSdt_gas_(i) = dSdt_gas_adv + dSdt_gas_src;

        if (pdata_->gas_.vr_(i) != 0.0 && pdata_->gas_.Sigma_(i) > pdata_->gas_.Sigma_floor_(i)) {
            // dt_temp = std::min(dt_temp, std::abs(pdata_->grid_.dr_(i) / pdata_->gas_.vr_(i)));
            dt_temp = std::min(dt_temp, pdata_->grid_.dr_(i) / (std::abs(pdata_->gas_.vr_(i)) + pdata_->gas_.cs_(i)));
        }

        if (dSdt_gas_src != 0.0 && pdata_->gas_.Sigma_(i) > pdata_->gas_.Sigma_floor_(i)) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->gas_.Sigma_(i)/dSdt_gas_src));
        }

        // dust
        dSdt_dust_adv = - 2.0*M_PI* (pdata_->grid_.r_bnd_(i+1)*flux_dust_(i+1) - pdata_->grid_.r_bnd_(i)*flux_dust_(i)) * dV_inv;
        dSdt_dust_src = pdata_->dust_.fdg_*(mdot_inf_r_(i+1) - mdot_inf_r_(i)) * dV_inv;
        dSdt_dust_(i) = dSdt_dust_adv + dSdt_dust_src;

        if (pdata_->dust_.vr_(i) != 0.0 && pdata_->dust_.Sigma_(i) != 0.0) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->grid_.dr_(i) / pdata_->dust_.vr_(i)));
        }

        if (dSdt_dust_src != 0.0 && pdata_->dust_.Sigma_(i) != 0.0) {
            dt_temp = std::min(dt_temp, std::abs(pdata_->dust_.Sigma_(i)/dSdt_dust_src));
        }

    }

    mdotloss_ = (1.0 + pdata_->dust_.fdg_)*mdot_inf_r_(is) - 2.0*M_PI*pdata_->grid_.r_bnd_(is)*(flux_gas_(is) + flux_dust_(is));

    return dt_temp;
}


void DiskEvolution::CalculateTimeStep(double dt)
{
    int i, j;
    int is = pdata_->grid_.is_, ie = pdata_->grid_.ie_;
    double sigmag_tot, sigmad_tot, dV;

    sigmag_tot = 0.0;
    sigmad_tot = 0.0;

    for (int i = is; i < ie; ++i) {

        pdata_->gas_.Sigma_(i) += dSdt_gas_(i) * dt;

        pdata_->dust_.Sigma_(i) += dSdt_dust_(i) * dt;

        dV = pdata_->grid_.r_vol_(i);
        sigmag_tot += pdata_->gas_.Sigma_(i) * dV;
        sigmad_tot += pdata_->dust_.Sigma_(i) * dV;
    }

    pdata_->star_.mass_ += mdotloss_ * dt;
    pdata_->mdot_ = mdotloss_;
    pdata_->gas_.total_mass_ = sigmag_tot;
    pdata_->dust_.total_mass_ = sigmad_tot;
    pdata_->total_disk_mass_ = sigmag_tot + sigmad_tot;
    pdata_->total_mass_ += ((1.0 + pdata_->dust_.fdg_)*mdot_inf_r_(ie) - mdot_wind_r_(ie)) * dt;
    pdata_->total_wind_loss_mass_ +=  mdot_wind_r_(ie)*dt;
    
    return;
}

