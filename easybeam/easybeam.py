"""
FEM Beam Analysis Tools
=========================

=================== ==========================================================
FEM Beam default tools
=============================================================================
Beam                Default 1D beam element
BeamSolver          Object that sticks beam elements together, and solve it
discertization      Optional discretization tool. By discretization you can
                    acquire precise results.
coordinates_array   Tool for creating array with coordinates of beams in tuple
moments_array       Tool for creating array with values of bending moments
                    of beams in tuple
shears_array        Tool for creating array with values of shears forces
                    of beams in tuple
disps_array         Tool for creating array with values of displacments
                    of beams in tuple
rotations_array     Tool for creating array with values of rotational angles
                    of beams in tuple
=================== ==========================================================

"""

import numpy as np
import numpy.linalg as lp
from itertools import chain

# My very first FEM solver program
# Author's name: Beniamin Dudek
# AGH UST WGiG
# field of studies: Civil Engineering


class No_Data(ValueError):
    pass


class Beam:
    # TODO:
    #   > use loads method \ad beginning
    #       to avoid bugs with discretization
    #   > define sample boundary at beginning
    # 1D Beam Finite Element

    def __init__(self, length, youngs_modulus=1, section=False):

        self.length = length  # unit: m
        self.youngs_modulus = youngs_modulus
        self.section = section  # units: m
        self.stifness()

        # some default definitons
        self.boundary(
            True,
            False,
            True,
            False
        )
        self.loads(
            linear_load=1000
        )

    def __repr__(self):

        return '{}({},{},{})'.format(
            __class__.__name__,
            self.length,
            self.youngs_modulus,
            self.section
            )

    # Stifness matrix for 1D beam element.
    def stifness(self):

        if self.section:
            self.area = self.section.area  # m^2
            self.moment_of_inertia_y = self.section.moment_of_inertia_y  # m^4
            self.moment_of_inertia_z = self.section.moment_of_inertia_z  # m^4
            self.elastic_modulus_y = self.section.elastic_modulus_y  # m^3
            self.elastic_modulus_z = self.section.elastic_modulus_z  # m^3
        else:
            self.area = 1  # unit: m^2
            self.moment_of_inertia_y = 1  # unit: m^4
            self.moment_of_inertia_z = 1  # unit: m^4
            self.elastic_modulus_y = 1  # unit: m^3
            self.elastic_modulus_z = 1  # unit m^3

        try:
            self.stifness_matrix = ((2*self.youngs_modulus*self.moment_of_inertia_z)/(self.length**3)) * np.array([[6, 3*self.length, -6, 3*self.length], [3*self.length, 2*self.length**2, -3*self.length, self.length**2], [-6, -3*self.length, 6, -3*self.length], [3*self.length, self.length**2, -3*self.length, 2*self.length**2]])
        except:
            raise No_Data

    def boundary(self, vertical_1, rotation_1, vertical_2, rotation_2):
        self.vertical_1 = bool(vertical_1)
        self.rotation_1 = bool(rotation_1)
        self.vertical_2 = bool(vertical_2)
        self.rotation_2 = bool(rotation_2)
        self.boundaries = [
            self.vertical_1,
            self.rotation_1,
            self.vertical_2,
            self.rotation_2
            ]

    def loads(self, force_1=0, force_2=0, moment_1=0, moment_2=0, linear_load=0):
        try:
            self.linear_load = linear_load

            # stworzone w celu ewentualnej dyskretyzacji
            self.init_force_1 = force_1
            self.init_force_2 = force_2
            self.init_moment_1 = moment_1
            self.init_moment_2 = moment_2

            self.force_1 = force_1 + self.linear_load*self.length*0.5
            self.force_2 = force_2 + self.linear_load*self.length*0.5
            self.moment_1 = moment_1 + self.linear_load*(self.length**2)*(1/12)
            self.moment_2 = moment_2 - self.linear_load*(self.length**2)*(1/12)

            self.system_loads = np.array(
                [self.force_1, self.moment_1, self.force_2, self.moment_2]
                )
        except:
            raise No_Data

    def internal_forces(self, solved_forces):
        self.__solved_forces__ = solved_forces
        self.internal_forces_array = \
            (-1)*np.transpose(np.matrix(self.system_loads)) + \
            self.stifness_matrix * np.transpose(
                np.matrix(self.__solved_forces__)
                )

        self.internal_forces_array = np.array(
            np.transpose(self.internal_forces_array)
            )[0]

    def solve(self):
        return self.__solved_forces__


class BeamSolver:

    def __init__(self, *beams):
        self.beams = tuple(chain.from_iterable(beams))

        self.n = 4 + (len(self.beams)-1)*2

        self.internal_agregation()
        self.internal_boundaries()
        self.internal_system_loads()
        self.apply_boundaries()
        self.solver()

    def internal_agregation(self):

        matrixes = [matrix.stifness_matrix for matrix in self.beams]

        self.global_stiffness_matrix = np.zeros([self.n, self.n])

        for i in range(len(matrixes)):
            current_matrix = np.zeros([self.n, self.n])
            current_matrix[(0+2*i):(4+2*i), (0+2*i):(4+2*i)] = matrixes[i]
            self.global_stiffness_matrix = \
                self.global_stiffness_matrix + current_matrix

    def internal_boundaries(self):

        vectors = [np.array(beam.boundaries) for beam in self.beams]

        self.global_boundary_vector = np.zeros([self.n])

        for i in range(len(vectors)):
            current_vector = np.zeros([self.n])
            current_vector[(0+2*i):(4+2*i)] = vectors[i]
            self.global_boundary_vector = \
                self.global_boundary_vector + current_vector

        self.global_boundary_vector = np.array(
            list(map(bool, self.global_boundary_vector))
            )

    def internal_system_loads(self):

        vectors = [beam.system_loads for beam in self.beams]

        self.global_system_loads_vector = np.zeros([self.n])

        for i in range(len(vectors)):
            current_vector = np.zeros([self.n])
            current_vector[(0+2*i):(4+2*i)] = vectors[i]
            self.global_system_loads_vector = \
                self.global_system_loads_vector + current_vector

    def apply_boundaries(self):

        self.boundariezed_stiffness_matrix = self.global_stiffness_matrix
        self.boundariezed_system_loads_vector = self.global_system_loads_vector
        self.__solved_vector__ = np.zeros([self.n])

        delete_counter = 0

        for i in range(np.size(self.global_boundary_vector, 0)):
            if self.global_boundary_vector[i]:
                self.boundariezed_stiffness_matrix = np.delete(
                    self.boundariezed_stiffness_matrix,
                    delete_counter,
                    0
                    )
                self.boundariezed_stiffness_matrix = np.delete(
                    self.boundariezed_stiffness_matrix,
                    delete_counter,
                    1
                    )
                self.boundariezed_system_loads_vector = np.delete(
                    self.boundariezed_system_loads_vector,
                    delete_counter
                    )
            else:
                delete_counter += 1
                self.__solved_vector__[i] = True

    def solver(self):
        self.unknowns_solved = lp.solve(
            self.boundariezed_stiffness_matrix,
            self.boundariezed_system_loads_vector
            )
        self.global_solvings_vector = np.zeros([self.n])

        unknowns_counter = 0
        for i in range(len(self.__solved_vector__)):
            if self.__solved_vector__[i]:
                self.global_solvings_vector[i] = \
                    self.unknowns_solved[unknowns_counter]
                unknowns_counter += 1

        for counter, beam_element in enumerate(self.beams):
            beam_element.internal_forces(
                self.global_solvings_vector[(0+2*counter):(4+2*counter)]
                )

    def results(self):

        return self.beams


def discretization(beam, number_of_elements):
    digitized_beam = []
    length_of_element = beam.length / number_of_elements

    for i in range(number_of_elements):
        # tworze tymczasowy element do dyskretyzacji
        # o klasie Beam
        current_beam = Beam(
            length_of_element,
            beam.youngs_modulus,
            beam.section
            )

        # okreslam warunki brzegowe poszczegolnych elementow
        # i nadaje obciazenia
        if i == 0:
            current_beam.boundary(beam.vertical_1, beam.rotation_1, False, False)
            current_beam.loads(
                force_1=beam.init_force_1,
                moment_1=beam.init_moment_1,
                linear_load=beam.linear_load
                )
        elif i == (number_of_elements-1):
            current_beam.boundary(
                False,
                False,
                beam.vertical_2,
                beam.rotation_2
                )
            current_beam.loads(
                force_2=beam.init_force_2,
                moment_2=beam.init_moment_2,
                linear_load=beam.linear_load
                )
        else:
            current_beam.loads(linear_load=beam.linear_load)
            current_beam.boundary(False, False, False, False)

        # dodaje gotowy element do listy
        digitized_beam.append(current_beam)

    return tuple(digitized_beam)


def momments_array(*beams):
    # fucntions for crating arrays with
    # values of bending moments
    beams = tuple(chain.from_iterable(beams))

    moments_array = []
    for beam in beams:
        moments = [
            beam.internal_forces_array[1],
            -1*beam.internal_forces_array[3]
            ]
        for moment in moments:
            moments_array.append(moment)

    return np.array(moments_array)


def shears_array(*beams):
    # fucntions for crating arrays with
    # values of shear forces
    beams = tuple(chain.from_iterable(beams))

    shears_array = []
    for beam in beams:
        shears = [
            beam.internal_forces_array[0],
            -1*beam.internal_forces_array[2]
            ]
        for shear in shears:
            shears_array.append(shear)

    return np.array(shears_array)


def disps_array(*beams):
    beams = tuple(chain.from_iterable(beams))

    disps_array = []
    for beam in beams:
        disps = [beam.solve()[0], beam.solve()[2]]
        for disp in disps:
            disps_array.append(disp)

    return np.array(disps_array)


def rotations_array(*beams):
    beams = tuple(chain.from_iterable(beams))

    rot_array = []
    for beam in beams:
        rots = [beam.solve()[1], beam.solve()[3]]
        for rot in rots:
            rot_array.append(rot)

    return np.array(rot_array)


def coordinates_array(*beams):
    # fucntions for crating arrays with
    # values of coordinates
    beams = tuple(chain.from_iterable(beams))
    coordinates_array = []

    current_length = 0
    for beam in beams:
        coordinates_array.append(current_length)
        current_length = current_length + beam.length
        coordinates_array.append(current_length)

    return np.array(coordinates_array)
