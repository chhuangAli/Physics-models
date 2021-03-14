class beam
{
    private:

        double time_duration = 0.;
        double beam_intensity = 0.;
        const double single_e_charge = 1.6e-19;
        const double mass_e = 9.1e-31;
        const double speed_light = 3.0e8;
        const double plank_const = 6.626e-34;

    public:

        double beam_energy = 0.;
        double beam_size[3] = {0}; // in x, y, z
        double number_of_particles_in_beam;
	double wavelength_gamma = 0.;
	double emittance_electron = 0.;

        void set_beam_duration(double);
        double get_beam_duration(void);
        void set_beam_size(double, double, double, double);
        double cal_num_elec_in_bunch(double);
        double cal_num_photon_in_pulse(double, double);
        double cal_beta_elec(double); // one input which is energy in eV.
	double cal_energy_laser(double);


};
