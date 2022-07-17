#pragma once
#include <cg/utils/ratio.h>
#include <cmath>

namespace cg {
#define Unit inline constexpr ratio

Unit dimensionless = 1.0;

/* Distance */
Unit f77unit = 1.0;
Unit angstrom = f77unit / 5.0;
Unit nanometer = 10.0 * angstrom;
Unit meter = nanometer * 1.0e9;

/* Time */
Unit nanosecond = 1.0;
Unit tau = 1.0 * nanosecond;
Unit microsecond = nanosecond * 1.0e3;
Unit millisecond = nanosecond * 1.0e6;
Unit second = nanosecond * 1.0e9;

/* Quantity */
Unit atom = 1.0;
Unit mol = 6.02214076e23 / atom;

/* Energy */
Unit eps = 1.0; /* \approx 1.5kcal/mol */
Unit kcal = eps * mol / 1.57824959;
Unit Joule = kcal / 4184.0;

/* Temperature */
Unit eps_kB = 1.0;
Unit kB = eps / eps_kB;
Unit Kelvin = 1.380649e-23 * Joule / kB;

/* Mass */
Unit kg = Joule * second * second / (meter * meter);
Unit amu = kg * 0.99999999965e-3 / mol;

/**
 * In the Fortran version of the model, distance of \p f77unit, time of
 * \p tau, energy of \p eps and the average mass of an aminoacid were units;
 * these are however incongruent, it's a confirmed bug.
 */
Unit f77mass = eps * tau * tau / (f77unit * f77unit);

/* EM stuff */
Unit echarge = 1.0;
Unit Coulomb = echarge / 1.602176634e-19;
Unit Ampere = Coulomb / second;
Unit cspeed = 299792458.0 * meter / second;
Unit Henry = kg * meter * meter / (second * second * Ampere * Ampere);
Unit mu_0 = 1.25663706212e-6 * Henry / meter;
Unit eps_0 = 1.0 / (mu_0 * cspeed * cspeed);

/* Degrees */
Unit rad = 1.0;
Unit deg = (2.0 * M_PI / 360.0) * rad;

#undef Unit
} // namespace cg