#! surface_file.py
#
# A class that represents a TICRATools Rotationally Symmetric Surface File
# for use in GRASP, CHAMP, etc.
#
import numpy
from scipy.interpolate import interp1d

class surface_of_rotation(object):
    """A class for representing and writing a TICRATools Rotationally Symmetric
    Surface File.

    Definition of the file format is given in the TICRATools Help PDF, p.2112.

    The rho array may be given in arbitrary order, with the i'th points of rho
    and z defining the i'th point.

    Attributes:
        header: A string for identification text, etc.
        kspace: A boolean indicating if the points are _not_ equispaced.
        np: An integer count of the number of points
        ktip: A boolean indicating if the surface is _not_ perpendicular to the
            rotation axis at the center.
        rho: 1-d numpy array holding the radii to points on the surface.
        z: 1-d numpy array holding the z coordinate of the points on the surface.
        rs: A float with the starting point of equispaced rho points.
        re: A float with the end point of the equispaced rho points.

    """
    def __init__(self, rho=None, z=None, rho_start=None, rho_end=None, np=None,
                auto_interpolate=False):
        header = "TICRA Tools Surface File"
        self.kspace = False
        self.ktip = False
        self._np = np
        self._rho = None
        self._z = None
        self._rs = rho_start
        self._re = rho_end
        self.auto_interpolate = auto_interpolate

        if rho != None:
            self.rho = rho
            self.kspace = True

        if z != None:
            self.z = z

        if z is None and np != None:
            # We have some number of points, so lets set up a placeholder z array
            # and fill in some sensible defaults for rs and re.
            self._np = np
            self._z = numpy.zeros((np), dtype=float)
            if self._rs is None:
                self._rs = 0.0
            if self._re is None:
                self._re = self._rs + 1.0

    @property
    def np(self):
        """Return the number of points in the z (and rho) arrays"""
        return self._np

    @np.setter
    def np(self, n):
        """Set the number of points in the surface.

        If the rho points are defined by the start and end points and a number
        of points, interpolate the z data to the new number of points."""
        if n is None:
            return
        if self.kspace:
            raise RuntimeWarning("Number of points is determined by the length of the z and rho arrays")
        else:
            self._np = n
            self._interpolate_z(self.rho)

    def is_consistent(self):
        """Checks surface_file object for self-consistency.

        Currently checks that a z array exists, and that if kspace is set, a
        rho array exists and is the same length as the z array.

        Returns: a boolean"""
        if self._z is None:
            return False

        if self.kspace:
            if self._rho is None:
                return False
            if len(self._rho) != len(self._z):
                return False

        return True

    @property
    def z(self):
        """Return the z array"""
        if not self.is_consistent():
            raise RuntimeWarning("surface_file object's z and rho arrays are inconsistent")
        return self._z

    @z.setter
    def z(self, z):
        """Set the surface z points.

        Also updates _np with number of points in z."""
        self._z = z
        self._np = len(z)

    @property
    def rho(self):
        """Returns the rho array.

        If points are evenly spaced, calculate the array from rs, re and np"""
        if self.kspace:
            return self._rho
        else:
            return numpy.linspace(self.rs, self.re, np)

    @rho.setter
    def rho(self, rho, interpolate=None):
        """Sets the rho array.  Also sets kspace to true if not already set"""
        if interpolate is None:
            interpolate = self.auto_interpolate

        self.kspace = True

        if interpolate:
            self._z = self._interpolate_z(rho=rho)

        self._rho = rho


    @property
    def rs(self):
        """Return the starting value of rho.

        If kspace is true, this is the smallest value in the rho array"""
        if self.kspace:
            return numpy.min(self._rho)
        else:
            return self._rs

    @rs.setter
    def rs(self, rs):
        """Set the starting value of rho."""
        self._rs = rs

    @property
    def re(self):
        """Return the ending value of rho.

        If kspace is true, this is the largest value in the rho array"""
        if self.kspace:
            return numpy.max(self._rho)
        else:
            return self._re

    @re.setter
    def re(self, re):
        """Set the ending value of rho"""
        self._re = re

    def _interpolate_z(self, rho):
        """Interpolate the z array to the points in the rho array"""
        z_interp = interp1d(self._rho, self._z)
        self._z = z_interp(rho, kind="cubic", fill_value="extrapolate")

    def write(self, filename):
        """Write out the surface to a TICRATools .rsf file"""
        if not self.is_consistent():
            raise RuntimeError("surface_file is not consistent!")

        f = open(filename, "w")

        f.write("{:s}\n".format(self.header))
        f.write("{:d}, {:d}, {:d}\n".format(self.np, self.kspace, self.ktip))
        if self.kspace:
            for i, r in enumerate(self.rho):
                f.write("{:f}, {:f}\n".format(r, self.z[i]))
        else:
            f.write("{:f}, {:f}\n".format(self.rs, self.re))
            for i, z in enumerate(self.z):
                f.write("{:f}\n".format(z))
        f.close()


# A convenience function to write a rho and z array to a given filename
def rho_z_to_file(rho, z, filename):
    """Write rho and z arrays to filename.

    Creates a surface_file opject"""
