# M-ANEOS

M-ANEOS is a FORTRAN program for the construction of thermodynamic equations of state, which extends the ANEOS computer code developed at Sandia National Laboratories [(Thompson and Lauson, 1972)](./docs/Thompson-Lauson-1972-Improvements-II.pdf).

## A brief description of ANEOS

ANEOS uses a suite of analytical approximations in different parts of thermodynamic phase space to construct an equation of state for use in shock physics calculations. The underlying analytical framework of ANEOS is an expression for the Helmholtz free energy in terms of its fundamental variables density and temperature, which is decomposed into three components: (a) a cold component, describing atomic and electronic interactions as a function of compression at absolute zero temperature; (b) a thermal component, describing the temperature-dependent parts of the interatomic forces, which vanishes at absolute zero temperature and approaches a perfect gas EOS at high pressures and low densities; and (c) an electronic component, describing ionization of electrons which becomes important at very high temperatures and low densities. This approach is advantageous because all useful thermodynamic functions can be derived from the Helmholtz free energy, and its derivatives, and similarly decomposed into separate components (cold, thermal, electronic). 

The full suite of analytical expressions used to define the various components of the thermodynamic functions in different parts of phase space are described in the ANEOS manual [(Thompson and Lauson, 1972)](./docs/Thompson-Lauson-1972-Improvements-II.pdf). The original version of ANEOS includes the functionality to account for three phase transitions: solid-liquid (melting); solid/liquid-vapor (vaporization); and a single solid-solid phase transition. However, the solid-liquid and solid-solid phase transitions could not be accounted for simultaneously.

## M-ANEOS: Modifications to ANEOS

[Melosh (2007)](./docs/Melosh2007.pdf) made significant improvements to ANEOS for the treatment of geologic materials. These improvements include the treatment of the vapor phase as composed of molecular clusters, rather than a monatomic mixture of atoms, and the use of a Mie-type potential in the expanded solid state, rather than a Morse potential. 

[Collins and Melosh (2014)](./docs/CollinsPosterLPSC2014.pdf) modified the ANEOS code to allow both a single solid-solid phase transition and the melt transition to be included simultaneously. 

[Stewart et al. (2019)](./docs/Stewart-2019-ANEOS_Modifications.pdf) modified ANEOS to include a parameter for the adjustment of the Debye thermal model used for the solid and liquid phases. The original ANEOS code had hard-coded a Dulong-Petit limit for the thermal term in the Helmholtz free energy [(Thompson and Lauson, 1972)](./docs/ANEOS_manual.pdf).

## Citing M-ANEOS and M-ANEOS-derived Equations of State in Scientific Works

To cite the M-ANEOS software we recommend including the following references:

Thompson SL, Lauson HS. [Improvements in the Chart D radiation-hydrodynamic CODE III: Revised analytic equation of state.](./docs/Thompson-Lauson-1972-Improvements-III.pdf) Albuquerque, N. Mex., USA: Sandia National Laboratories; 1972.

Melosh HJ. [A hydrocode equation of state for SiO2](./docs/Melosh2007.pdf). Meteoritics & Planetary Science 2007; 42:2079–98.

Thompson SL, Lauson HS, Melosh HJ, Collins GS and Stewart, ST, M-ANEOS: A Semi-Analytical Equation of State Code, Zenodo, http://doi.org/10.5281/zenodo.3525030

To cite a specific M-ANEOS-derived equation of state, we recommend citing the source of the M-ANEOS input parameters, together with the three references listed above, using a format such as:

“In this work we used the forsterite equation of state (Stewart et al., 2019) constructed using M-ANEOS (Thompson and Lauson, 1972; Melosh, 2007; Thompson et al., 2019).”

An exemplar for the documentation of an M-ANEOS derived equation of state is:

Stewart, Sarah T., Davies, Erik J., Duncan, Megan S., Lock, Simon J., Root, Seth, Townsend, Joshua P., Kraus, Richard G., Caracas, Razvan, Jacobsen, Stein B., 2019. Equation of State Model Forsterite-ANEOS-SLVTv1.0G1: Documentation and Comparisons. [https://doi.org/10.5281/zenodo.3478631]

## M-ANEOS Contributors

S. L. Thompson

H. S. Lauson

H. J. Melosh

G. S. Collins

S. T. Stewart

## Installation and usage

M-ANEOS is provided with a very bare-bones makefile, which should be modified for your fortran compiler.

To compile, run `make` in the `src` directory. This creates library `libaneos.a` and executable `m-aneos`.

The command `make install` will install `m-aneos` in the example directory, in addition to an example M-ANEOS input file for quartz (see below).

To run the example, change in to the `example/` directory and invoke `./m-aneos`. This will use the example `ANEOS.INPUT` file and the provided `tablegrid.txt` file to produce SESAME tables in two different formats, as well as the standard ANEOS output file `ANEOS.OUTPUT`.

To run M-ANEOS without generating an output table, invoke `./m-aneos --no_table`. This will produce the standard ANEOS output file `ANEOS.OUTPUT` only.

Simple usage instructions are printed with `./m-aneos --help`.

Available input files can be found in the `input/` directory.

## Input files

quartz.input - Melosh H. J. (2007) [A hydrocode equation of state for SiO2](./docs/Melosh2007.pdf). Meteoritics & Planetary Science, 42:2079–98.

serpentine.input - Brookshaw, L. (1998) [An Equation of State for Serpentine](./docs/Brookshaw.pdf), Tech. Rep. Working Paper Series SC-MC-9813 (Queensland: Faculty of Sciences, University of Southern Queensland)

dunite.input - Collins G. S. and Melosh H. J. (2014) "[Improvements to ANEOS for multiple phase transitions](./docs/CollinsPosterLPSC2014.pdf)," in 45th Lunar and Planetary Science Conference, p. [2664](./docs/CollinsANEOSLPSC2014.pdf).
