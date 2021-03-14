#include "beam.h"
#include <cmath>

void beam::set_beam_duration(double time)
{
    time_duration = time;
}

double beam::get_beam_duration(void)
{
    return time_duration;
}


void beam::set_beam_size(double x_size, double y_size, double time, double speed)
{
    beam_size[0] = x_size;
    beam_size[1] = y_size;
    time_duration = time;
    beam_size[2] = (time_duration * speed); // for z size
//    return longitudinal_size;
}

double beam::cal_num_elec_in_bunch(double intensity )
{
    double number_of_elec_in_bunch = 0.;

    //take care of unit of the intensity
    beam_intensity = intensity;

    number_of_elec_in_bunch = beam_intensity / single_e_charge;
    return number_of_elec_in_bunch;
}

double beam::cal_num_photon_in_pulse(double intensity, double wave_length)
{
    double num_photons_in_laser = 0.;

    beam_intensity = intensity;
    double E_photon = 0.;

    E_photon = plank_const * speed_light / wave_length;
    num_photons_in_laser = beam_intensity / E_photon;

    return num_photons_in_laser;
}

double beam::cal_beta_elec(double energy)
{
    double gamma = 0.;
    beam_energy = energy*single_e_charge; // change to joules
    gamma = (beam_energy)/(mass_e*speed_light*speed_light);
    double beta = 0.; // output
    beta = sqrt(1.-(1./pow(gamma,2)));

    return beta;
}

double beam::cal_energy_laser(double wavelength)
{
  wavelength_gamma = wavelength;
  beam_energy = plank_const * speed_light / wavelength_gamma;

  return beam_energy;
}
