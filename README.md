# planet_heat
Simulation of a purely rocky planet's heat, with an aim to eventually include parameters such as rotation rate, axial tilt, seasons, and elliptical orbits, probably in that order.

This first simulation had the aim of purely simulating the surface of the planet near to a star using the differential equation

Laplacian(r,theta) T(R,theta) = constants * T(R,theta)^4 .

This equation is non-linear so normal separation of variables doesn't work, but I hoped that using spherical coordinates and only looking at the temperature 
on the surface at a radius R, I could simply solve this equation in one dimension, then fit those points to the interior with Bessel functions and spherical 
harmonics. This did not work though. The first simualtion is kept here for posterity.
