# residence-time-acf

* This program calculates residence time autocorrelation functions to evaluate exchange kinetics 
in and out of a defined region.

* This program is inspired from the methodology described in:
S. A. Kislenko, R. H. Amirov, I. S. Samoylov, Phys. Chem. Chem. Phys., 12, 11245 (2010) 

* The program calculates the autocorrelation functions for two heaviside functions, hI and hC:

Intermittent function = <hI(t).hI(0)>/<hI(0)^2>

Continuous function = <hC(t).hC(0)>/<hC(0)^2>

* To compile the program you can use any C compiler, e.g.:

  `gcc -o resid-time-acf.x resid-time-acf.c -lm`

* To execute the program, you just have to launch it with an input file using:

  `./resid-time-acf.x < inputfile`
  
  The input file should provide the following information

1. Name of the xyz file with the positions 
   
   At this point the program only works for a constant number of atoms ordered in the same way 
in all configurations.

2. Number of configurations in the positions file

3. Sampling parameter for configurations, i.e. configurations are considered every nskip
   If nskip = 1, all configurations are considered.

4. Time between two considered configurations

   Note that no conversion is done so the time unit in the output is the same as the input.

5. Number of species for which the calculation will be done

   Note that there is currently a limit of 1000 atom types.

6. Atom type for each species (one line per species)

7. Limits in the z direction 

8. Number of correlations steps
