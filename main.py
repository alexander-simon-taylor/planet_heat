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

        
    def tidal_locked_temp(self, Bessel_num, sph_harm_num, min_steps):
        if min_steps < (Bessel_num + 1) * (sph_harm_num + 1):
            print("Too few steps.")
            return 0
        else:
            if np.mod(min_steps, 2) == 0:
                min_steps += 1
            temperature_vector = np.zeros(((Bessel_num + 1) * (sph_harm_num + 1), 1))
            for i in range(1, (Bessel_num + 1) * (sph_harm_num + 1)):
                temperature_vector[i] = np.pi * i / ((Bessel_num + 1) * (sph_harm_num + 1) - 1)
            print(temperature_vector)
            RK4_array = np.zeros((min_steps, 3))
            RK4_array[0] = [initial_temperature, initial_gradient, 0]
            temp_vec_counter = 0
            step_size = np.pi / (min_steps - 1)
            for i in range(1, min_steps):
                RK4_array[i] = self.RK4_step(RK4_array[i - 1, 0], RK4_array[i - 1, 1], step_size, i * step_size)
                if i * step_size >= temperature_vector[temp_vec_counter]:
                    temperature_vector[temp_vec_counter] = RK4_array[i, 0]
                    temp_vec_counter += 1
            print(temperature_vector)
            print(RK4_array)
            return 0

    def tidal_lock_B_gradient(self, temperature, temp_gradient, theta):
        B_gradient = stefan_boltzmann * self.planet_heat_diffusivity * self.planet_radius ** 2 * temperature ** 4
        print(B_gradient)
        B_gradient /= self.planet_density * self.planet_shc
        B_gradient += temp_gradient / np.tan(theta)
        if theta < np.pi / 2:  # only adds the star's heat if the side of the planet is facing the star
            B_gradient_addon = self.planet_heat_diffusivity * self.planet_radius ** 2 * self.star_luminosity * np.cos(
                theta)
            B_gradient_addon /= 4 * np.pi * self.star_distance ** 2 * self.planet_density * self.planet_shc
        else:
            print("Gets here")
            B_gradient_addon = 0
        print(B_gradient)
        print(B_gradient_addon)
        return B_gradient_addon - B_gradient

    def RK4_step(self, temperature, temp_gradient, step_size, theta):
        # q is the temperature coordinate, k is the gradient coordinate
        q1 = step_size * temp_gradient
        k1 = step_size * self.tidal_lock_B_gradient(temperature, temp_gradient, theta)

        q2 = temp_gradient + k1 / 2
        k2 = step_size * self.tidal_lock_B_gradient(temperature + q1 / 2, q2, theta)
        q2 *= step_size

        q3 = temp_gradient + k2 / 2
        k3 = step_size * self.tidal_lock_B_gradient(temperature + q2 / 2, q3, theta)
        q3 *= step_size

        q4 = temp_gradient + k3
        k4 = step_size * self.tidal_lock_B_gradient(temperature + q3, q4, theta)
        q4 *= step_size

        new_temperature = temperature + (q1 + 2 * q2 + 2 * q3 + q4) / 6
        new_gradient = temp_gradient + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return [new_temperature, new_gradient, theta]


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
