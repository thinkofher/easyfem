# Simple Sections Classes
# for Fem_Beam Numerical Calculations
# AUTHOR: Beniamin Dudek

from math import pi
from copy import deepcopy


class Section():

    def __init__(self, moment_of_inertia_y, moment_of_inertia_z, area):

        self.moment_of_inertia_y = moment_of_inertia_y
        self.moment_of_inertia_z = moment_of_inertia_z
        self.area = area

    def return_strong_i_moment(self):

        return self.moment_of_inertia_y

    def return_weak_i_moment(self):

        return self.moment_of_inertia_z

    def return_area(self):

        return self.area


class RectangleSection(Section):

    def __init__(self, height, width):

        self.width = width
        self.height = height

        if self.width == 0 or self.height == 0:
            raise(ValueError)

        self.area = self.height*self.width
        self.moments_of_inertias()
        self.elastic_moduluses()

    def __repr__(self):

        return '{}({}, {})'.format(__class__.__name__, self.height, self.width)

    def moments_of_inertias(self):

        self.moment_of_inertia_y = (self.height**3 * self.width) / 12
        self.moment_of_inertia_z = (self.height * self.width**3) / 12

    def elastic_moduluses(self):

        self.elastic_modulus_y = self.moment_of_inertia_y * (self.height/2)
        self.elastic_modulus_z = self.moment_of_inertia_z * (self.width/2)

    def rotate(self):

        self.__initial_values__ = ( deepcopy(self.height), deepcopy(self.width) )
        self.height = self.__initial_values__[1]
        self.width = self.__initial_values__[0]

        self.moments_of_inertias()
        self.elastic_moduluses


class CircleSection(Section):

    def __init__(self, diameter):

        self.diameter = diameter

        if self.diameter == 0:
            raise(ValueError)

        self.area = pi * (self.diameter/2)**2
        self.moments_of_inertias()
        self.elastic_moduluses()

    def __repr__(self):

        return '{}({})'.format(__class__.__name__, self.diameter)

    def moments_of_inertias(self):

        self.moment_of_inertia_y = pi * (self.diameter**4 / 64)
        self.moment_of_inertia_z = self.moment_of_inertia_y

    def elastic_moduluses(self):

        self.elastic_modulus_y = self.moment_of_inertia_y * (self.diameter/2)
        self.elastic_modulus_z = self.moment_of_inertia_z * (self.diameter/2)


class HollowCircle(Section):

    def __init__(self, big_diamater, small_diamater):

        self.big_diamater = big_diamater
        self.small_diamater = small_diamater
        self.outside_circle = CircleSection(self.big_diamater)
        self.inside_circle = CircleSection(self.small_diamater)
        self.area = self.outside_circle.area - self.inside_circle.area

        self.moments_of_inertias()
        self.elastic_moduluses()

    def __repr__(self):

        return '{}({},{})'.format(
            __class__.__name__,
            self.big_diamater,
            self.small_diamater
            )

    def moments_of_inertias(self):

        self.moment_of_inertia_y = self.outside_circle.moment_of_inertia_y - \
            self.inside_circle.moment_of_inertia_y
        self.moment_of_inertia_z = self.outside_circle.moment_of_inertia_z - \
            self.inside_circle.moment_of_inertia_z

    def elastic_moduluses(self):

        self.elastic_modulus_y = self.outside_circle.elastic_modulus_y - \
            self.inside_circle.elastic_modulus_y
        self.elastic_modulus_z = self.outside_circle.elastic_modulus_z - \
            self.inside_circle.elastic_modulus_z


class HollowRectangle(Section):

    def __init__(self, height, width, thickness):

        self.outside_height = height
        self.outside_width = width
        self.thickness = thickness
        self.inside_height = height - 0.5*self.thickness
        self.inside_width = width - 0.5*self.thickness

        self.outside_rectangle = RectangleSection(
            self.outside_height,
            self.outside_width
            )
        self.inside_rectangle = RectangleSection(
            self.inside_height,
            self.inside_width
            )
        self.area = self.outside_rectangle.area - self.inside_rectangle.area

        self.moments_of_inertias()
        self.elastic_moduluses()

    def __repr__(self):

        return '{}({},{},{})'.format(
            __class__.__name__,
            self.outside_height,
            self.outside_width,
            self.thickness
            )

    def moments_of_inertias(self):

        self.moment_of_inertia_y = \
            self.outside_rectangle.moment_of_inertia_y - \
            self.inside_rectangle.moment_of_inertia_y
        self.moment_of_inertia_z = \
            self.outside_rectangle.moment_of_inertia_z - \
            self.inside_rectangle.moment_of_inertia_z

    def elastic_moduluses(self):

        self.elastic_modulus_y = \
            self.outside_rectangle.elastic_modulus_y - \
            self.outside_rectangle.elastic_modulus_y
        self.elastic_modulus_z = \
            self.outside_rectangle.elastic_modulus_z - \
            self.outside_rectangle.elastic_modulus_z


class IBeamSection(Section):

    def __init__(self, height, width, flange_thickness, web_thickness):

        self.height = height
        self.width = width
        self.flange_thickness = flange_thickness
        self.web_thickness = web_thickness

        self.area = 2 * self.flange_thickness * self.width + \
            (self.height - 2 * self.flange_thickness) * self.web_thickness
        self.moments_of_inertias()
        self.elastic_moduluses()

    def __repr__(self):

        return '{}({}, {}, {}, {})'.format(
                                            __class__.__name__,
                                            self.height,
                                            self.width,
                                            self.flange_thickness,
                                            self.web_thickness
                                            )
    def moments_of_inertias(self):

        self.moment_of_inertia_y = (1/12) * ( self.width * self.height**3 - ( self.width - self.web_thickness ) * ( self.height - 2 * self.flange_thickness )**3 )
        self.moment_of_inertia_z = (1/12) * ( 2 * self.flange_thickness * self.width**3 + ( self.height - 2 * self.flange_thickness )*self.web_thickness**3 ) 

    def elastic_moduluses(self):

        self.elastic_modulus_y = self.moment_of_inertia_y * (self.height/2)
        self.elastic_modulus_z = self.moment_of_inertia_z * (self.width/2)
