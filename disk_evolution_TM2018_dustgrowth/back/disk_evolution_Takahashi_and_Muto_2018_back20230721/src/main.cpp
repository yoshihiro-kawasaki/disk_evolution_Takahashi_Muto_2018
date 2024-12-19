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


int main(int argc, char* argv[])
{
    std::string input_file;
    if (argc != 2) {
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

    // output time
    // double t       = 0.0;                     
    // double tend    = 7.7e5*pc::SOLAR_YEAR;
    // double dt, dtmin, dt_disk, dt_collpse;
    // Array1D<double> time(7);
    // time(0) = 1.5e5 * pc::SOLAR_YEAR;
    // time(1) = 1.8e5 * pc::SOLAR_YEAR;
    // time(2) = 2.65e5 * pc::SOLAR_YEAR;
    // time(3) = 3.5e5 * pc::SOLAR_YEAR;
    // time(4) = 4.3e5 * pc::SOLAR_YEAR;
    // time(5) = 6.0e5 * pc::SOLAR_YEAR;
    // time(6) = tend;

    double t       = 0.0;                     
    double tend    = 2.0e5*pc::SOLAR_YEAR;
    double dt, dtmin, dt_disk, dt_collpse;
    Array1D<double> time(5);
    time(0) = 1.2e5 * pc::SOLAR_YEAR;
    time(1) = 1.4e5 * pc::SOLAR_YEAR;
    time(2) = 1.6e5 * pc::SOLAR_YEAR;
    time(3) = 1.8e5 * pc::SOLAR_YEAR;
    time(4) = 2.0e5 * pc::SOLAR_YEAR;

    dt      = 1e80;
    dtmin   = 1e80;
    dt_disk = 1e80;

    double mdot_inf, r_ini,rho_ini,drinidt,t_ini, c_wind;
    mdot_inf  = 0.0;
    drinidt   = 0.0;

    // set initial condition
    r_ini = be_sphere.CalculateInitialRinit();
    if (r_ini >= be_sphere.GetRoutBonner()) {
        mdot_inf = 0.0;
        drinidt  = 0.0;
        c_wind = sim.gas_.c_wind_;
    } else {
        be_sphere.CalculateMdotInfall(r_ini, drinidt, mdot_inf);
        c_wind = 0.0;
    }
    t_ini = sqrt(0.5*CUB(r_ini)/(pc::GRAVITATIONAL_CONSTANT*sim.star_.mass_init_))*2.56717; // 2.56717は圧力勾配力の影響
    t = t_ini;
    gevol.SetInitialCondtion();
    gevol.CalculateDiskGasParamters(mdot_inf, r_ini, c_wind);

    int count = 0;
    int time_index = 0;
    std::string outfilename;
    bool switch_infall = false;

    // main integration
    while (t < tend) {

        count += 1;

        // determine dt
        gevol.CalculateFlux();
        if (drinidt == 0.0) {
            dt_collpse = 1.0e2 * pc::ASTRONOMICAL_UNIT;
        } else {
            dt_collpse = 1.0e-2*r_ini / drinidt;
        }
        dt_disk = 3.0e-1*gevol.CalculateDtMin();
        dtmin = std::min(dt_collpse, dt_disk);
        dtmin = std::min(3.0e0*pc::SOLAR_YEAR, dtmin);

        if ((time(time_index) - t) <= dtmin) {
            dt = time(time_index) - t;
        } else {
            dt = dtmin;
        }
        // if (t == t_ini) dt = 1.0e-10 * pc::SOLAR_YEAR;

        if (count % 50000 == 0) {
            std::cout << count << " " << (dt/pc::SOLAR_YEAR) << " " << (t/pc::SOLAR_YEAR) << std::endl;
        }

        // time step (n → n + 1)
        gevol.CalculateTimeStep(dt);
        r_ini  += drinidt*dt;
        t += dt;
        sim.time_ = t;

        if (r_ini >= be_sphere.GetRoutBonner()) {
            mdot_inf = 0.0;
            drinidt  = 0.0;
            c_wind = sim.gas_.c_wind_;
            if (switch_infall == false) {
                std::cout << "end infall : time = " << (t/pc::SOLAR_YEAR) << std::endl;
                switch_infall = true;
            }
        } else {
            be_sphere.CalculateMdotInfall(r_ini, drinidt, mdot_inf);
            c_wind = 0.0;
        }
        gevol.CalculateDiskGasParamters(mdot_inf, r_ini, c_wind);

        // output
        if (t >= time(time_index)) {
            outfilename = sim.outdir_ + "/disk" + std::to_string(time_index);
            sim.Output(outfilename);
            time_index += 1;
            std::cout << std::scientific << "output : t = " << (t/pc::SOLAR_YEAR) << " solar year" << std::endl;
        }

    }


    std::cout << std::scientific;
    std::cout << "star mass            = " << (sim.star_.mass_/pc::SOLAR_MASS) << " solar mass" << std::endl; 
    std::cout << "total disk gas mass  = " << (sim.gas_.total_mass_/pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total disk dust mass = " << (sim.dust_.total_mass_/pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total disk mass      = " << (sim.total_disk_mass_/pc::SOLAR_MASS) << " solar mass"  << std::endl;
    std::cout << "total mass           = " << (sim.total_mass_/pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "wind mass loss       = " << (sim.total_wind_loss_mass_/pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "m                    = " << ((sim.total_mass_ + sim.total_wind_loss_mass_)/pc::SOLAR_MASS) << " solar mass" << std::endl;
    std::cout << "total count          = " << count << std::endl;


    end = std::chrono::system_clock::now();
    const double timed = static_cast<double>(std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count());
    std::cout << "time = " << (timed*1.0e-3) << "sec" << std::endl;
    

    return 0;
}