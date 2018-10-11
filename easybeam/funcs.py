"""
easyfem Analysis Tools
=========================

=================== ==========================================================
easyfem funcs tools
=============================================================================
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
from itertools import chain
from easyfem.easybeam.classes import Beam

# My very first FEM solver program
# Author's name: Beniamin Dudek
# AGH UST WGiG
# field of studies: Civil Engineering


def discretization(beam, number_of_elements):
    digitized_beam = []
    length_of_element = beam.length / number_of_elements

    for i in range(number_of_elements):
        # creating teporary element for discretization
        current_beam = Beam(
            length_of_element,
            beam.youngs_modulus,
            beam.section
            )

        # declaring boundaries and loads for each element
        if i == 0:
            current_beam.boundary(
                beam.vertical_1,
                beam.rotation_1,
                False, False
                )
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

        # adding ready element to list
        digitized_beam.append(current_beam)

    return tuple(digitized_beam)


def momments_array(*beams):
    '''
    fucntions for crating arrays with
    values of bending moments
    '''
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
    '''
    fucntions for crating arrays with
    values of shear forces
    '''
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
    '''
    fucntions for crating arrays with
    values of displacemnets
    '''
    beams = tuple(chain.from_iterable(beams))

    disps_array = []
    for beam in beams:
        disps = [beam.solve()[0], beam.solve()[2]]
        for disp in disps:
            disps_array.append(disp)

    return np.array(disps_array)


def rotations_array(*beams):
    '''
    fucntions for crating arrays with
    values of rotations
    '''
    beams = tuple(chain.from_iterable(beams))

    rot_array = []
    for beam in beams:
        rots = [beam.solve()[1], beam.solve()[3]]
        for rot in rots:
            rot_array.append(rot)

    return np.array(rot_array)


def coordinates_array(*beams):
    '''
    fucntions for crating arrays with
    values of coordinates
    '''
    beams = tuple(chain.from_iterable(beams))
    coordinates_array = []

    current_length = 0
    for beam in beams:
        coordinates_array.append(current_length)
        current_length = current_length + beam.length
        coordinates_array.append(current_length)

    return np.array(coordinates_array)
