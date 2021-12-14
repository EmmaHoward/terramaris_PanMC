# (C) British Crown Copyright 2010 - 2017, Met Office
#
# This file is part of Iris.
#
# Iris is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Iris is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with Iris.  If not, see <http://www.gnu.org/licenses/>.
"""
Calculus operations on :class:`iris.cube.Cube` instances.

See also: :mod:`NumPy <numpy>`.

"""

#from __future__ import (absolute_import, division, print_function)
from six.moves import (filter, input, map, range, zip)  # noqa
import six

import re

import cf_units
import numpy as np

from iris._deprecation import warn_deprecated
import iris.cube
import iris.coords
import iris.coord_systems
import iris.analysis
import iris.analysis.maths
from iris.analysis.cartography import (DEFAULT_SPHERICAL_EARTH_RADIUS,
                                       DEFAULT_SPHERICAL_EARTH_RADIUS_UNIT)
from iris.util import delta


__all__ = ['cube_delta', 'differentiate', 'curl']


def delta(ndarray, dimension, circular=False):
    """
    Calculates the difference between values along a given dimension.

    Args:

    * ndarray:
        The array over which to do the difference.

    * dimension:
        The dimension over which to do the difference on ndarray.

    * circular:
        If not False then return n results in the requested dimension
        with the delta between the last and first element included in
        the result otherwise the result will be of length n-1 (where n
        is the length of ndarray in the given dimension's direction)

        If circular is numeric then the value of circular will be added
        to the last element of the given dimension if the last element
        is negative, otherwise the value of circular will be subtracted
        from the last element.

        The example below illustrates the process::

            original array              -180, -90,  0,    90
            delta (with circular=360):    90,  90, 90, -270+360

    .. note::

        The difference algorithm implemented is forward difference:

            >>> import numpy as np
            >>> import iris.util
            >>> original = np.array([-180, -90, 0, 90])
            >>> iris.util.delta(original, 0)
            array([90, 90, 90])
            >>> iris.util.delta(original, 0, circular=360)
            array([90, 90, 90, 90])

    """
    if circular is not False:
        _delta1 = np.roll(ndarray, -1, axis=dimension)
        _delta2 = np.roll(ndarray, 1, axis=dimension)
        last_element = [slice(None, None)] * ndarray.ndim
        last_element[dimension] = slice(-1, None)
        first_element = [slice(None, None)] * ndarray.ndim
        first_element[dimension] = slice(None, 1)
#
        if not isinstance(circular, bool):
            result = np.where(ndarray[last_element] >= _delta1[last_element])[0]
            _delta1[last_element] -= circular
            _delta1[last_element][result] += 2*circular
            result = np.where(ndarray[first_element] >= _delta2[first_element])[0]
            _delta2[first_element] -= circular
            _delta2[first_element][result] += 2*circular
#
#
        np.subtract(_delta1, _delta2, _delta1)
    else:
        forward = [slice(None,None)] * ndarray.ndim
        backward = [slice(None,None)] * ndarray.ndim
        forward[dimension] = slice(2,None)
        backward[dimension] = slice(None,-2)
        _delta1 = np.subtract(ndarray[forward],ndarray[backward])
#
    return _delta1/2



def _construct_delta_coord(coord):
    """
    Return a coordinate of deltas between the given coordinate's points.
    If the original coordinate has length n and is circular then the
    result will be a coordinate of length n, otherwise the result will be
    of length n-2.

    """
    if coord.ndim != 1:
        raise iris.exceptions.CoordinateMultiDimError(coord)
    circular = getattr(coord, 'circular', False)
    if coord.shape == (1,) and not circular:
        raise ValueError('Cannot take interval differences of a single '
                         'valued coordinate.')
    if coord.shape == (2,) and not circular:
        raise ValueError('Cannot take centred interval differences of a two '
                         'valued coordinate.')
    if circular:
        circular_kwd = coord.units.modulus or True
    else:
        circular_kwd = False

    if coord.bounds is not None:
        bounds = delta(coord.bounds, 0, circular=circular_kwd)
    else:
        bounds = None

    points = delta(coord.points, 0, circular=circular_kwd)
    new_coord = iris.coords.AuxCoord.from_coord(coord).copy(points, bounds)
    new_coord.rename('change_in_%s' % new_coord.name())

    return new_coord


def _construct_midpoint_coord(coord, circular=None):
    """
    Return a coordinate of mid-points from the given coordinate. If the
    given coordinate has length n and the circular flag set then the
    result will be a coordinate of length n, otherwise the result will be
    of length n-1.

    """
    if circular and not hasattr(coord, 'circular'):
        raise ValueError('Cannot produce circular midpoint from a coord '
                         'without the circular attribute')

    if circular is None:
        circular = getattr(coord, 'circular', False)
    elif circular != getattr(coord, 'circular', False):
        warn_deprecated('circular flag and Coord.circular attribute do '
                        'not match')

    if coord.ndim != 1:
        raise iris.exceptions.CoordinateMultiDimError(coord)
    if coord.shape == (1,) and not circular:
        raise ValueError('Cannot take the midpoints of a single valued '
                         'coordinate.')

    # Calculate the delta of the coordinate
    # (this deals with circularity nicely).
    mid_point_coord = _construct_delta_coord(coord)

    # if the coord is circular then include the last one, else, just take 0:-1
    if circular:
      circular_slice = slice(0, None)
    else:
      circular_slice = slice(1, -1)
#    import pdb;pdb.set_trace()
    if coord.bounds is not None:
#        axis_delta = mid_point_coord.bounds
        mid_point_bounds = coord.bounds[circular_slice, :]
    else:
        mid_point_bounds = None

    # Get the deltas
    axis_delta = mid_point_coord.points
#    assert (axis_delta[0] == axis_delta).all()
    # Add half of the deltas to the original points
    # if the coord is circular then include the last one, else, just take 0:-1
    mid_point_points = coord.points[circular_slice]

    # Try creating a coordinate of the same type as before, otherwise,
    # make an AuxCoord.
    try:
        mid_point_coord = coord.from_coord(coord).copy(mid_point_points,
                                                       mid_point_bounds)
    except ValueError:
        mid_point_coord = iris.coords.AuxCoord.from_coord(coord).copy(
            mid_point_points, mid_point_bounds)

    return mid_point_coord


def cube_delta(cube, coord):
    """
    Given a cube calculate the difference between each value in the
    given coord's direction.


    Args:

    * coord
        either a Coord instance or the unique name of a coordinate in the cube.
        If a Coord instance is provided, it does not necessarily have to
        exist in the cube.

    Example usage::

        change_in_temperature_wrt_pressure = \
cube_delta(temperature_cube, 'pressure')

    .. note:: Missing data support not yet implemented.

    """
    # handle the case where a user passes a coordinate name
    if isinstance(coord, six.string_types):
        coord = cube.coord(coord)

    if coord.ndim != 1:
        raise iris.exceptions.CoordinateMultiDimError(coord)

    # Try and get a coord dim
    delta_dims = cube.coord_dims(coord.name())
    if ((coord.shape[0] == 1 and not getattr(coord, 'circular', False)) or
            not delta_dims):
        raise ValueError('Cannot calculate delta over {!r} as it has '
                         'length of 1.'.format(coord.name()))
    delta_dim = delta_dims[0]

    # Calculate the actual delta, taking into account whether the given
    # coordinate is circular.
    delta_cube_data = delta(cube.data, delta_dim,
                            circular=getattr(coord, 'circular', False))

    # If the coord/dim is circular there is no change in cube shape
    if getattr(coord, 'circular', False):
        delta_cube = cube.copy(data=delta_cube_data)
    else:
        # Subset the cube to the appropriate new shape by knocking off
        # the last row of the delta dimension.
        subset_slice = [slice(None, None)] * cube.ndim
        subset_slice[delta_dim] = slice(1, -1)
        delta_cube = cube[tuple(subset_slice)]
        delta_cube.data = delta_cube_data

    # Replace the delta_dim coords with midpoints
    # (no shape change if circular).
    for cube_coord in cube.coords(dimensions=delta_dim):
        delta_cube.replace_coord(_construct_midpoint_coord(
            cube_coord, circular=getattr(coord, 'circular', False)))

    delta_cube.rename('change_in_{}_wrt_{}'.format(delta_cube.name(),
                                                   coord.name()))

    return delta_cube


def differentiate(cube, coord_to_differentiate):
    r"""
    Calculate the differential of a given cube with respect to the
    coord_to_differentiate.

    Args:

    * coord_to_differentiate:
        Either a Coord instance or the unique name of a coordinate which
        exists in the cube.
        If a Coord instance is provided, it does not necessarily have to
        exist on the cube.

    Example usage::

        u_wind_acceleration = differentiate(u_wind_cube, 'forecast_time')

    The algorithm used is equivalent to:

    .. math::

        d_i = \frac{v_{i+1}-v_i}{c_{i+1}-c_i}

    Where ``d`` is the differential, ``v`` is the data value, ``c`` is
    the coordinate value and ``i`` is the index in the differential
    direction. Hence, in a normal situation if a cube has a shape
    (x: n; y: m) differentiating with respect to x will result in a cube
    of shape (x: n-1; y: m) and differentiating with respect to y will
    result in (x: n; y: m-1). If the coordinate to differentiate is
    :attr:`circular <iris.coords.DimCoord.circular>` then the resultant
    shape will be the same as the input cube.

    In the returned cube the `coord_to_differentiate` object is
    redefined such that the output coordinate values are set to the
    averages of the original coordinate values (i.e. the mid-points).
    Similarly, the output lower bounds values are set to the averages of
    the original lower bounds values and the output upper bounds values
    are set to the averages of the original upper bounds values. In more
    formal terms:

    * `C[i] = (c[i] + c[i+1]) / 2`
    * `B[i, 0] = (b[i, 0] + b[i+1, 0]) / 2`
    * `B[i, 1] = (b[i, 1] + b[i+1, 1]) / 2`

    where `c` and `b` represent the input coordinate values and bounds,
    and `C` and `B` the output coordinate values and bounds.

    .. note:: Difference method used is the same as :func:`cube_delta`
        and therefore has the same limitations.

    .. note:: Spherical differentiation does not occur in this routine.

    """
    # Get the delta cube in the required differential direction.
    # This operation results in a copy of the original cube.
    delta_cube = cube_delta(cube, coord_to_differentiate)

    if isinstance(coord_to_differentiate, six.string_types):
        coord = cube.coord(coord_to_differentiate)
    else:
        coord = coord_to_differentiate

    delta_coord = _construct_delta_coord(coord)
    delta_dim = cube.coord_dims(coord.name())[0]

    # calculate delta_cube / delta_coord to give the differential.
    delta_cube = iris.analysis.maths.divide(delta_cube, delta_coord, delta_dim)

    # Update the standard name
    delta_cube.rename('derivative_of_{}_wrt_{}'.format(cube.name(),
                                                       coord.name()))
    return delta_cube

