#include "array.hpp"
#include "physical_constants.hpp"
#include "sim_data.hpp"
#include "disk_evolution.hpp"
#include "bonner_ebert_sphere.hpp"

#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <sys/stat.h>

int main(int argc, char *argv[])
{
    std::string input_file;
    if (argc != 2)
    {
        std::cerr << "# ERROR : main, no input file" << std::endl;
        std::exit(1);
    }

    input_file = argv[1];
    std::cout << "input file : " << input_file << std::endl;

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();

    Utils::InputConfigure input;
    std::string input_file_type;
    input.ReadFile(input_file);

    //
    SimulationData sim(input);
    sim.Initialization();
    // return 0;

    BonnerEbertSphere be_sphere(&sim);
    be_sphere.SetBonnerEbertSphere();
    // be_sphere.Output("be.txt");

    DiskEvolution gevol(&sim);

    double t = 0.0;
    double tend = 2.0e5 * pc::SOLAR_YEAR;
    double dt, dtmin, dt_disk, dt_collpse, tout;
    Array1D<double> time(300);
    for (int i = 0; i < 300; ++i) {
        time(i) = (1.0 + 0.1 * double(i)) * 1.0e5 * pc::SOLAR_YEAR;
    }

    dt = 1e80;
    dtmin = 1e80;
    dt_disk = 1e80;

    double mdot_inf, r_ini, rho_ini, drinidt, t_ini, c_wind;
    mdot_inf = 0.0;
    drinidt = 0.0;

    // set initial condition
    r_ini = be_sphere.CalculateInitialRinit();
    if (r_ini >= be_sphere.GetRoutBonner())
    {
        mdot_inf = 0.0;
        drinidt = 0.0;
        c_wind = sim.gas_.c_wind_;
    }
    else
    {
        be_sphere.CalculateMdotInfall(r_ini, drinidt, mdot_inf);
        c_wind = 0.0;
    }
    t_ini = sqrt(0.5 * CUB(r_ini) / (pc::GRAVITATIONAL_CONSTANT * sim.star_.mass_init_)) * 2.56717; // 2.56717は圧力勾配力の影響
    t = t_ini;
    tout = t + 5.0e1 * pc::SOLAR_YEAR;
    gevol.SetInitialCondtion();
    gevol.CalculateDiskGasParamters(mdot_inf, r_ini, c_wind);

    int count = 0;
    int time_index = 0;
    std::string outfilename;
    bool switch_infall = false;

    std::string logfile = sim.outdir_ + "/log";
    std::ofstream file(logfile, std::ios::out | std::ios::trunc);
    FAILED_TO_OPEN(file, logfile);

    // main integration
    while (t < tend)
    {

        count += 1;

        // determine dt
        gevol.CalculateFlux();
        if (drinidt == 0.0) {
            dt_collpse = 1.0e2 * pc::ASTRONOMICAL_UNIT;
        } else {
            dt_collpse = 1.0e-2 * r_ini / drinidt;
        }
        dt_disk = 3.0e-1 * gevol.CalculateDtMin();
        dtmin = std::min(dt_collpse, dt_disk);
        if (dtmin == 0.0) {
            std::cout << dt_disk << " " << dt_collpse << std::endl;
        }
        dtmin = std::min(3.0e0 * pc::SOLAR_YEAR, dtmin);

        if ((time(time_index) - t) <= dtmin) {
            dt = time(time_index) - t;
        } else {
            dt = dtmin;
        }

        if (count % 200000 == 0) {
            std::cout << count << " " << (dt / pc::SOLAR_YEAR) << " " << (t / pc::SOLAR_YEAR) << std::endl;
        }

        // time step (n → n + 1)
        gevol.CalculateTimeStep(dt);
        r_ini += drinidt * dt;
        t += dt;
        sim.time_ = t;

        if (r_ini >= be_sphere.GetRoutBonner()) {
            mdot_inf = 0.0;
            drinidt = 0.0;
            c_wind = sim.gas_.c_wind_;
            if (switch_infall == false) {
                std::cout << "end infall : time = " << (t / pc::SOLAR_YEAR) << std::endl;
                switch_infall = true;

                // outfilename = sim.outdir_ + "/disk_end_infall";
                outfilename = sim.outdir_ + "/disk" + std::to_string(count);
                sim.OutputData(outfilename);
                sim.LogOutput(count, file);
                // break;
            }
        } else {
            be_sphere.CalculateMdotInfall(r_ini, drinidt, mdot_inf);
            c_wind = 0.0; //sim.gas_.c_wind_;
        }
        gevol.CalculateDiskGasParamters(mdot_inf, r_ini, c_wind);

        // output
        if (t >= tout || t >= time(time_index)) {

            // outfilename = sim.outdir_ + "/disk" + std::to_string(time_index);
            outfilename = sim.outdir_ + "/disk" + std::to_string(count);
            // sim.Output(outfilename);
            sim.OutputData(outfilename);
            sim.LogOutput(count, file);

            if (t >= time(time_index)) {
                time_index += 1;
                std::cout << std::scientific << "output : t = " << (t / pc::SOLAR_YEAR) << " solar year" << std::endl;
            } else {
                if (switch_infall == false) {
                    tout += 5e1 * pc::SOLAR_YEAR;
                } else {
                    tout += 2.0e2 * pc::SOLAR_YEAR;
                }
            }
        }
    }

    std::cout << std::scientific;
    std::cout << "star mass            = " << (sim.star_.mass_ / pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total disk gas mass  = " << (sim.gas_.total_mass_ / pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total disk dust mass = " << (sim.dust_.total_mass_ / pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total disk mass      = " << (sim.total_disk_mass_ / pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total mass           = " << (sim.total_mass_ / pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total infall mass    = " << (sim.total_infall_mass_ / pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "wind mass loss       = " << (sim.total_wind_loss_mass_ / pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "m                    = " << ((sim.total_mass_ + sim.total_wind_loss_mass_) / pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total count          = " << count << std::endl;

    file.close();

    end = std::chrono::system_clock::now();
    const double timed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "time = " << (timed * 1.0e-3) << "sec" << std::endl;

    return 0;
}