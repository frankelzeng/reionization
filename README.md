# Fall 2018 Reionization

This is the main repository for Chenxiao's fall 2018 research project focusing on post-reionization temperature.

The code takes the ionization front velocity and reionization source temperature as inputs, and calculate the post-reionization temperaturefor five species (e, HI, HII, HeI, HeII) in the grid model. This new calculation include the non-equilibrium interaction bewteen species and study how the post-I temperature is shifted from the equilibrium case when different reionization parameters are applied

Use "cc -std=c99 -Wall -o gif gif_190127_2.c" for compiling

For the front velocity 5e8 cm/s and source blackbody temperature 5e4 K, we run the executable with: "./gif 50000 500000000 > output.txt". The numbers can be replaced by other desired input parameters. 

The output file will generate a table of post-reionization temperature for each species with 5 columns and 2000 rows (number of grids in the model).

Rerun the calculation with gif_180817_original.c to figure out the post-I temperature for the equilibrium case

The energy transfer rate between species are interpolated from experiments and the detail calculations are shown in the mathematica notebook reionization_cgs.nb
