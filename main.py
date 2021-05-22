import numpy as np
# import math

stefan_boltzmann = 5.67e-8
initial_temperature = 395
initial_gradient = 0


class planet_star_system:
    star_luminosity = 3.84e26;
    star_distance = 149000000000;
    planet_shc = 830;
    planet_radius = 6371000;
    planet_density = 2650;
    planet_rotation_rate = 0;
    planet_axial_tilt = 0;
    planet_orbital_rate = 0;
    planet_heat_diffusivity = 0;

    def __init__(self, init_star_luminosity, init_star_distance, init_planet_shc,
                 init_planet_radius, init_planet_density, init_rotation_rate,
                 init_axial_tilt, init_orbital_rate, init_planet_diffusivity):
        self.star_luminosity = init_star_luminosity
        self.star_distance = init_star_distance
        self.planet_shc = init_planet_shc
        self.planet_radius = init_planet_radius
        self.planet_density = init_planet_density
        self.planet_rotation_rate = init_rotation_rate
        self.planet_axial_tilt = init_axial_tilt
        self.planet_orbital_rate = init_orbital_rate
        self.planet_heat_diffusivity = init_planet_diffusivity

    def create_grid(self, grid_size, import_static = None):
        # This creates a grid of temperature values. A cube grid is used to make indexing easy, despite this program
        # only considering a spherical region.
        if grid_size//2 == 0:
            grid_size += 1
            # Makes grid size an odd number, to ensure a centre point
        grid_radius = (grid_size + 1)/2 + 1
        # The +1 at the end forces there to be at least one grid point of empty space around the planet sphere
        if import_static:
            # To-do, add ability to unpack list of Fourier coefficients to create starting grid
        else:
            temperature_grid = np.zeros((grid_size + 2,grid_size + 2,grid_size + 2))
            for x in range(-grid_radius,grid_radius): # x is axis that faces towards the star
                for y in range(-grid_radius,grid_radius):
                    previous_empty_space = True
                    for z in range(-grid_radius,grid_radius):
                        if x**2 + y**2 + z**2 > grid_radius**2: # If point is in empty space
                            temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius] = np.array([0,"space"])
                            previous_empty_space = True
                        else: # If point is the surface or inside the surface
                            if previous_empty_space:
                                temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius] = np.array([0, "surface"])
                            else:
                                temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius] = np.array([0, "interior"])
                            previous_empty_space = False

star_luminosity = 3.84e26;
star_distance = 1.49e11;
planet_shc = 830;
planet_radius = 6371000;
planet_density = 2650;
planet_rotation_rate = 0;
planet_axial_tilt = 0;
planet_orbital_rate = 0;
planet_heat_diffusivity = 1e-6;

Earth_Sun = planet_star_system(star_luminosity, star_distance, planet_shc,
                               planet_radius, planet_density, planet_rotation_rate,
                               planet_axial_tilt, planet_orbital_rate, planet_heat_diffusivity)
Earth_Sun.tidal_locked_temp(3, 3, 35)