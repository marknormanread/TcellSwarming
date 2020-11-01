package utils;

import java.util.Arrays;

import sim.util.Double2D;
import sim.util.Double3D;

/**
 * This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
 * @author Mark N. Read
 *
 */
public class Utils 
{
	public static double radiansToDegrees(double rad)
	{	return rad * 180.0 / Math.PI; 	}


	public static double degrees_to_radians(double deg)
	{    return deg * Math.PI / 180.0;	}
	
	/** Computes the dot product between two vectors, represented as Double3D objects. 
	 */
	public static double dotProduct(Double3D a, Double3D b)
	{
		return (a.x*b.x + a.y*b.y + a.z*b.z);
	}
	
	public static double dotProduct(Double2D a, Double2D b)
	{
		return (a.x*b.x + a.y*b.y);
	}
	
	/** Computes the cross product of two vectors. The cross product is a vector perpendicular to the plane in which 
	 *  vectors u and v lie. Its magnitude is equivalent to the area of the parallelogram with u and v as its edges.
	 *  The direction in which the cross product extends (there are two ways it can extend out of the plane) depends
	 *  on which orientation (right/left hand rule) describes the geometry of space, and the order in which vectors
	 *  are applied. i.e., cross product is non-commutative. 
	 *  
	 *  Returns a vector of zero magnitude if the two inputs are parallel. You may need to check for this. 
	 */
	public static Double3D crossProduct(Double3D u, Double3D v)
	{
		double nx = u.y*v.z - u.z*v.y;
		double ny = u.z*v.x - u.x*v.z;
		double nz = u.x*v.y - u.y*v.x;
		
		return new Double3D(nx, ny, nz);		
	}
	
	/**
	 * Returns a positive number of a, b, and c adhere to the right-hand rule or orientation/standard basis.
	 * Returns a negative number otherwise.
	 * Returns zero when all three vectors lie on the same plane.
	 */
	public static double tripleProduct(Double3D a, Double3D b, Double3D c)
	{
		Double3D cp = crossProduct(b, c);
		double dp = dotProduct(a, cp);
		return dp;
	}
		
	/** Projects vector 'w' ONTO vector 'v'. Neither w nor v need to be unit
	 *  vectors, the code below works for vectors of any length. 
	 */	
	public static Double3D projectVector(Double3D w, Double3D v)
	{
		double wDotv = dotProduct(w, v);
		double vDotv = dotProduct(v, v);
		double scalar = wDotv / vDotv;
		
		return new Double3D(v.x * scalar, v.y * scalar, v.z * scalar);
	}
	
	/** Calculates the angle between two vectors, in RADIANS. This is done given the rule:
	 *  a.b = |a| . |b| . cos(angle)
	 *  ie, dot product of a and b equals their lengths multiplied with one another and the cos of the angle between
	 *  the vectors. 
	 */
	public static double angleBetweenVectors(Double3D a, Double3D b)
	{
		double dp = dotProduct(a, b);
		double cosAng = dp / (length(a) * length(b));
		// this can happen because of precision errors. Math.acos can't take abs(value) > 1.0 
		// usually occurs when the vectors are parallel but pointing in opposite directions.  
		if (Math.abs(cosAng) > 1.0)	 
			cosAng /= cosAng;		// this preserves magnitude (+/-) and scales to 1.0;
		return Math.acos(cosAng);	// return value is in radians. 
	}
	
	/**
	 * The anticlockwise rotation around n required to get from vector a to vector b.
	 * n: the normal of the plane through which the rotation takes place.
	 * radians: if false, will return in degrees
	 * shortest: False, angle will be an anticlockwise rotation (can be > pi radians). If True, angles > pi radians
	 * are instead returned as negative values.
	 */
	public static double rotation_anticlockwise(Double3D a, Double3D b, Double3D n, boolean radians, boolean shortest)
	{
		double angle = angleBetweenVectors(a, b);
	
	    // if positive, then the three vectors are oriented according to the standard basis functions i, j, k.
	    double tp = tripleProduct(a, b, n);
	    if (tp < 0)
	        if (shortest)
	            angle = -angle;
	        else
	            angle = (Math.PI * 2.) - angle;

	    return radians? angle : radiansToDegrees(angle);
	}
	
	
	/** Calculates the length of vector 'w', often written |w|.  */
	public static double length(Double3D w)
	{
		return Math.sqrt(w.x*w.x + w.y*w.y + w.z*w.z);
	}
	public static double lengthSq(Double3D w)
	{
		return w.x*w.x + w.y*w.y + w.z*w.z;
	}
	
	/** Calculates the length of the vector described by the three arguments, 
	 *  which are distances along the three orthogonal axes.
	 */
	public static double length(double dx, double dy, double dz)
	{
		return Math.sqrt( dx*dx + dy*dy + dz*dz);
	}
	
	/** Calculates the displacement between the two points specified by the locations w and v. */
	public static double displacement(Double3D w, Double3D v)
	{
		return length(w.x - v.x, w.y - v.y, w.z - v.z);
	}
	
	/** Returns the vector that takes you FROM location w TO location v. */
	public static Double3D difference(Double3D w, Double3D v)
	{	return new Double3D(v.x - w.x, v.y - w.y, v.z - w.z);	}
	
	/** Returns the unit vector of 'w'. If 'w' has zero length, returns a zero vector.  */
	public static Double3D unitVector(Double3D w)
	{		
		double wLength = length(w);
		if (wLength != 0.0)
			return new Double3D(w.x/wLength, w.y/wLength, w.z/wLength);
		return w;		// return zero vector, which w already is. 
	}	
	
	/**
	 * Calculates the median value from an array. Array does not need to be sorted in advance. 
	 * 
	 * This code was taken from http://stackoverflow.com/questions/8938235/java-sort-an-array on 19/05/2014.
	 * It has been slightly modified by Mark Read. 
	 */
	public static double median(double[] m) 
	{
		Arrays.sort(m);				// algorithm requires array to be sorted. 
	    int middle = m.length/2;	// this floors down to the lower int. 
	    if (m.length%2 == 1) {
	        return m[middle];
	    } else {
	        return (m[middle-1] + m[middle]) / 2.0;
	    }
	}
}
