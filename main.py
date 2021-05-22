import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
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
♠
    def create_grid(self, grid_size, import_static = None):
        # This creates a grid of temperature values. A cube grid is used to make indexing easy, despite this program
        # only considering a spherical region.
        if grid_size//2 == 0:
            grid_size += 1
            # Makes grid size an odd number, to ensure a centre point
        grid_radius = int((grid_size + 1)/2)
        # The +1 at the end forces there to be at least one grid point of empty space around the planet sphere
        if import_static:
            return "Finished here" # To-do, add ability to unpack list of Fourier coefficients to create starting grid
        else:
            temperature_grid = np.zeros((grid_size + 2, grid_size + 2, grid_size + 2, 2))
            # First loop to categorise points into either space or interior points
            for x in range(-grid_radius, grid_radius):
                for y in range(-grid_radius, grid_radius):
                    for z in range(-grid_radius, grid_radius): # z is the axis facing towards the star
                        if x**2 + y**2 + z**2 > grid_radius**2: # If point is in empty space
                            temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius] = [0, 0] # 0 for space
                        else: # If point is the surface or inside the surface
                            temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius] = [0, 1] # 1 for interior
            # Second loop to find surface points, those being interior points next to at least one space point,
            # and give them an approximate temperature
            number_of_surface_points = 0
            for x in range(-grid_radius,grid_radius):
                for y in range(-grid_radius,grid_radius):
                    for z in range(-grid_radius,grid_radius): # z is the axis facing towards the star
                        if x**2 + y**2 + z**2 <= grid_radius**2 and x**2 + y**2 + z**2 > (grid_radius - 2)**2:
                            # If point is in rough boundary of surface
                            surface_check = False
                            if temperature_grid[x + grid_radius + 1, y + grid_radius, z + grid_radius, 1] == 0:
                                surface_check = True
                            elif temperature_grid[x + grid_radius - 1, y + grid_radius, z + grid_radius, 1] == 0:
                                surface_check = True
                            elif temperature_grid[x + grid_radius, y + grid_radius + 1, z + grid_radius, 1] == 0:
                                surface_check = True
                            elif temperature_grid[x + grid_radius, y + grid_radius - 1, z + grid_radius, 1] == 0:
                                surface_check = True
                            elif temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius + 1, 1] == 0:
                                surface_check = True
                            elif temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius - 1, 1] == 0:
                                surface_check = True
                            if surface_check:
                                number_of_surface_points += 1
                                if x >= 0:
                                    temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius - 1] = [CMB_temp, 2] # 2 for surface
                                else:
                                    starting_surface_temp = np.power(star_luminosity * np.abs(x) / (4 * np.pi * stefan_boltzmann * star_distance ** 2 * grid_radius) + CMB_temp ** 4, 1/4)
                                    temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius - 1] = [starting_surface_temp, 2]
            print(number_of_surface_points)
        return temperature_grid

    def surface_plot(self, temperature_grid):
        grid_radius = int((len(temperature_grid[:, 0, 0, 0]) - 1)/2)
        coordinate_points = np.zeros((0, 4))
        for x in range(-grid_radius, grid_radius):
            for y in range(-grid_radius, grid_radius):
                for z in range(-grid_radius, grid_radius):  # z is the axis facing towards the star
                    if temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius, 1] == 2:
                        print(z)
                        coordinate_row = [x, y, z, temperature_grid[x + grid_radius, y + grid_radius, z + grid_radius, 0]]
                        coordinate_points = np.vstack((coordinate_points, coordinate_row))
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        print(coordinate_points[:, 3])
        ax.scatter(coordinate_points[:, 0], coordinate_points[:, 1], coordinate_points[:, 2], marker=".", c=coordinate_points[:, 3])
        plt.show()
        return

star_luminosity = 3.84e26;
star_distance = 1.49e11;
planet_shc = 830;
planet_radius = 6371000;
planet_density = 2650;
planet_rotation_rate = 0;
planet_axial_tilt = 0;
planet_orbital_rate = 0;
planet_heat_diffusivity = 1e-6;
CMB_temp = 2.7;

Earth_Sun = planet_star_system(star_luminosity, star_distance, planet_shc,
                               planet_radius, planet_density, planet_rotation_rate,
                               planet_axial_tilt, planet_orbital_rate, planet_heat_diffusivity)
first_grid = Earth_Sun.create_grid(35)
Earth_Sun.surface_plot(first_grid)
