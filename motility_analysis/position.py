"""
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Mark N. Read, 2016.


Captures a single spatio-temporal observation. These can be collected together into tracks.

Mark N. Read, 2017
"""


class Position:
    """
    Position of a cell in time and space. Values can be assigned None if the track was temporarily lost. For
    instance, IMARIS can do this, it will recognise that two segments belong to the same track (because of temporal
    and spatial proximity) but doesn't have any data for the period the track was lost. In the case of IMARIS,
    this can happen if the cell's signal is particularly weak at that time.
    """
    def __init__(self, x, y, z, time_s, timeID):
        """
        :param x: Position on x axis. None if this is a place-holder object for an agent that was temporarily lost
        :param y: Position on y axis
        :param z: Position on z axis
        :param time_s: time stamp (in absolute real time) of this spatio-temporal position in seconds
        :param timeID: ID of temporal observation (ie, not absolute time, but integer designating observation
        number, starting from zero at the beginning of the experiment.
        """
        # these values are drawn straight from the input data
        self.x = x
        self.y = y
        self.z = z
        # the timestamp (in absolute real time) of this spatio-temporal position in seconds
        self.time_s = float(time_s)
        self.time_m = self.time_s / 60.
        self.timeID = timeID
        # these values have to be post-processed
        self.displacement = None  # how far the cell moved from previous time step to this one.
        self.displacementVector = None  # tuple (x,y,z) of displacement in each axis from previous step to this one
        self.instant_fmi = None  # instantaneous forward migration index.
        self.total_displacement = 0.0  # displacement from track's starting point to the position at this time.
        self.total_displacement_squared = 0.0  # displacement from track's starting point squared.
        self.turn = None  # the angular velocity through which the cell pitched from prev time step to this one.
        self.roll = None  # the angular velocity through which the cell rolled from prev time step to this one.
        self.speed = None  # speed of cell, from last time step to this one.
        # used to filter out data corresponding to cell "shudder" (cell essentially stopped).
        self.meets_arrest_coeff_threshold = True
        self.speed_acceleration = None
        self.turn_acceleration = None
        # set to True if this position is part of an interpolation to fix broken tracks.
        self.interpolated = False

    def tracked(self):
        """
        Returns true if the cell's location at this time was known. Position objects are created even for times when
        the cell position is not known. However, this has implications on some statistics, hence this method.
        """
        return self.x is not None and self.y is not None and self.z is not None

    def __str__(self):
        if self.x is not None:
            s = 'Position: x = {:.2f}, y = {:.2f}, z = {:.2f}, time = {:.2f} seconds'\
                .format(self.x, self.y, self.z, self.time_s)
        else:
            s = 'Position: x, y, = None (unknown), time = {:.2f} seconds'\
                .format(self.time_s)
        if self.speed is not None:
            s = s + ', speed = {:.2f}'.format(self.speed)
        else:
            s = s + ', speed = None'
        return s
