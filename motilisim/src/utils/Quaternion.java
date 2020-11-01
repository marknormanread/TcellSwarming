package utils;

import core.Simulation;
import sim.util.Double3D;


/**
 * 
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
 * 
 * Quaternions are seriously confusing. Read this before proceeding. Remember, we're talking about rotations of the 
 * x-axis vector (1,0,0) here, not vectors describing heading. 
 * 
 * Some notes on using Quaternions. Currently it is used to maintain the orientation of cells in space, relative to 
 * absolute coordinate system. Cell's are assumed to face their x-axes in relative space, and hence always move forward
 * along their x-axis. Rotations around their x-axis correspond to the cell rolling. Rotations around their y-axes
 * correspond to changes in their pitch.  
 * 
 * If the rotation of a cell relative to the absolute coordinate system is stored in 'orientation', then 
 * 'orientation.transform([1.0,0.0,0.0])' will give the vector in absolute space along which the cell currently faces. 
 * 
 * Successive rotations are applied to orientation to represent the cell's changing orientation in space. THE ORDER IN 
 * WHICH ROTATIONS ARE APPLIED MATTERS. Orientation changes applied relative to the cell, such as rolling it, or 
 * changing its pitch, must be applied thus:
 * 
 *   Quaternion change = Quaternion.representRotation(angle, 1, 0, 0) // roll cell along its x-axis, direction it faces.
 *   orientation = orientation.multiply(change). 
 * 
 * Hence, 'change' is added atop 'orientation'. 
 * 
 * The alternative is to specify a rotation in absolute space. For instance, when a cell has collided with something, 
 * and must bounce off it. This happens in this simulation, when cells collide with one another or with hair follicles. 
 * A vector in absolute space, 'bounce', represents normal of the contact. Here the rotation must be calculated in
 * absolute space (because bounce is represented in absolute space). 
 * 
 *   facing = orienation.transform([1,0,0])  // retrieve vector of direction cell faces, in absolute space. 
 * 	 normal = Utils.crossProduct(facing,bounce) // vector is the normal to the plane along which facing and bounce lie.
 * 					// rotation must be applied around this axis, such that facing, bounce and facing-after all lie
 * 					// in the same plane. 
 *   angle = xxx // this can be the angle between bounce and facing, or several other things, depending on whether the
 *   				// cell shoud slide along the obstacle, or bounce off it. 
 *   Quaternion bounceRot = Quaternion.representRotation(angle, [normal.x, normal.y, normal.z])
 *   orientation = bounceRot.mutliply(orientation);
 *   
 * In this example the order of multiplications is changed. Since the rotation and (importantly!) the rotation axis are
 * expressed in absolute space, they go first. 
 * 
 * Hence, the order in which the muliplicants are multiplied depends on whether the rotation is expressed in relative
 * or absolute space. 
 * 
 * TO MAKE A CELL FACE A PARTICULAR DIRECTION:
 * specify the rotation required in relation to the x-axis (along which cells face), and set that as the cell's 
 * orientation. Hence, you need to calculate the angle between the desired heading and the x axis, and express that 
 * rotation around the normal of the two vectors. 
 *  
 * Some simple examples, which do not consider the roll of the cell, which is also changed. 
 * Make the cell face the x-axis:
 *   orientation = Quaternion.representRotation(0., 1., 0., 0.);
 * 
 * Make the cell face along the negative x axis:
 *   orientation = Quaternion.representRotation(Math.PI, 0.0, 1.0, 0.0);
 *   (take the x axis (1,0,0) and apply a 180 degree rotation around the y axis. This could be done for z axis too, and
 *   the roll orientation of the cell would be different. )
 *   
 * Make the cell face along the y axis:
 *   orientation = Quaternion.representRotation(Math.PI / 2.0, 0.0, 0.0, 1.0); 
 *   (rotate 90 degrees around z axis.)
 * 
 * Make the cell face along z axis:
 *   orientation = Quaternion.representRotation(Math.PI/2, 0.0, 1.0, 0.0);
 *   
 * 
 * Code in this class is based on a number of very helpful tutorials, both accessed on 09/07/2013:
 * - http://www.cprogramming.com/tutorial/3d/quaternions.html
 * - http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/transforms/index.htm * 
 * 
 * 
 * @author Mark Read
 */
public class Quaternion 
{
	private double w;	// this is the rotation. 
	private double x;	// these are the vector components around which rotation occurs. 
	private double y;
	private double z;
	// used in checking when the Quaternion needs to be normalised. Errors in Quaternion maths compound.  
	private final double normalisation_tolerance = 0.00001;
	
	/** Constructs a Quaternion with the specified scalar and complex components. Do not use this constructor to
	 *  directly translate from conventional 3D space into 4D Quaternion space, there is another static method for that
	 *  (i.e., simply inputing an angle and a 3D vector to this constructor will not yield correct results in 
	 *  computations, the 3D representation needs to be transformed into 4D quaternion space.)
	 * 
	 * @param w	Scalar value, representing a rotation 
	 * @param x Complex components.  
	 * @param y
	 * @param z
	 */
	public Quaternion(double w, double x, double y, double z)
	{
		this.w = w;
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public static Quaternion identity()
	{
		return new Quaternion(1., 0., 0., 0.);
	}
	
	/** The magnitude of this Quaternion
	 * 
	 * @return The magnitude of this Quaternion.
	 */
	public double length()
	{	
		return Math.sqrt(w*w + x*x + y*y + z*z);
	}
	
	/** Returns the magnitude of this Quaternion squared. 
	 * 
	 * @return The magnitude of this Quaternion squared. 
	 */
	public double length2()
	{	
		return w*w + x*x + y*y + z*z;	
	}
	
	/** Quanternion maths (when implemented) can accumulate compounded errors as calculations are made. This performs
	 *  normalisation if it is needed. It modifies this Quanternion. 
	 * 
	 * @return This Quanternion object, now normalised. 
	 */
	public Quaternion normalise()
	{
		double len = this.length();
		w = w / len;
		x = x / len;
		y = y / len;
		z = z / len;		
		return this;
	}
	
	/** Returns the conjugate of this Quanternion. This is similar to complex numbers, where the conjugate is simply
	 *  the same Quanternion with its complex part (x y and z components) negated. A new Quanternion is returned to 
	 *  preserve this object.  
	 */
	public Quaternion conjugate()
	{
		return new Quaternion(this.w, -this.x, -this.y, -this.z);
	}
	
	/** Returns the inverse of this Quaternion as a new object, the original is preserved. The inverse can be used to 
	 *  undo a transformation. 
	 *  
	 *  b = q.transform(a)
	 *  q.inverse().transform(b) = a. 
	 */
	public Quaternion inverse()
	{	
		double len2 = this.length2();
		return new Quaternion(this.w/len2, -this.x/len2, -this.y/len2, -this.z/len2);
	}
	
	/** Multiplies this quaternion by that quaternion. Note that quaternion multiplication is non-commutative, hence
	 *  this*that does not equal that*this. This operation returns a new Quaternion, hence both this and that remain 
	 *  intact.  
	 *  
	 *  Compound rotations through quaternion mulitplication. If Rz = R1 followed by R2, to derive Rz you must calculate
	 *  both R1 and R2 as quaternions, and then multiple R2 by R1 (ie, reverse order). 
	 *  Rz = R2.multiple(R1). 
	 *  
	 * @param that The Quaternion to be multiplied as this*that.  
	 * @return The result of this*that.
	 */
	public Quaternion multiply(Quaternion that)
	{
		double new_w = this.w*that.w - this.x*that.x - this.y*that.y - this.z*that.z; 
		double new_x = this.w*that.x + this.x*that.w + this.y*that.z - this.z*that.y;
		double new_y = this.w*that.y - this.x*that.z + this.y*that.w + this.z*that.x;
		double new_z = this.w*that.z + this.x*that.y - this.y*that.x + this.z*that.w;
		
		return new Quaternion(new_w, new_x, new_y, new_z);
	}
	
	/** Converts the specified angle of rotation around the supplied axis (ie, in angle-axis form) into a Quanternion.
	 * 
	 * @param angle Angle of rotation, in radians. 
	 * @param x X component of vector around which rotation is applied.
	 * @param y Y component of vector around which rotation is applied.
	 * @param z Z component of vector around which rotation is applied.
	 * @return A new Quanternion representing the rotation in Quanternion space. 
	 */
	public static Quaternion representRotation(double angle, double x, double y, double z)
	{
		if (x*x + y*y + z*z != 1.0)
		{
			double len = Math.sqrt(x*x + y*y + z*z);
			x /= len;
			y /= len;
			z /= len;
		}
		double a = angle / 2.;
		double new_w = Math.cos(a);
		double new_x = x * Math.sin(a);
		double new_y = y * Math.sin(a);
		double new_z = z * Math.sin(a);
		return new Quaternion(new_w, new_x, new_y, new_z);
	}
	
	/** Returns a quaternion randomly selected from a uniform distribution in 4D (rotation) space. 
	 *  This is based on http://planning.cs.uiuc.edu/node198.html, accessed on 14/07/2013.  
	 */
	public static Quaternion randomUniform()
	{
		double u1 = Simulation.instance.random.nextDouble();
		double u2 = Simulation.instance.random.nextDouble();
		double u3 = Simulation.instance.random.nextDouble();
		Quaternion q = new Quaternion(
			Math.sqrt(1.0-u1) * Math.sin(2.0 * Math.PI * u2), 	// w.
			Math.sqrt(1.0-u1) * Math.cos(2.0 * Math.PI * u2), 	// x.
			Math.sqrt(u1) * Math.sin(2.0 * Math.PI * u3),		// y.
			Math.sqrt(u1) * Math.cos(2.0 * Math.PI * u3));		// z.
		return q;
	}	
	
	/** Returns a quaternion that if set as a cell's orientation will make the cell face the desired vector. 
	 * Note that this will only work if you use the x-axis as the direction a cell ordinarily faces. 
	 * If this is not the case, then you need to change the reference (currently xaxis) vector accordingly. 
	 * 
	 * @param dir direction in which the cell shuold face. 
	 */
	public static Quaternion faceVector(Double3D dir)
	{
		final Double3D xaxis = new Double3D(1.0, 0.0, 0.0);  
		// the normal is perpendicular to the plane on which bounce and facing vectors lie. 
		Double3D normal = Utils.crossProduct(xaxis, dir);
		if (normal.lengthSq() == 0.0)	// happens when dir is parallel to xaxis. 
			normal = new Double3D(0.0,1.0,0.0);	// select y axis as default, but any vector perpendicular to x will do.
		double angle = Utils.angleBetweenVectors(xaxis, dir);	// in radians.				
		return Quaternion.representRotation(angle, normal.x, normal.y, normal.z);
	}

	/** Applies the transformation represented by this Quaternion to the specified vector. For example, this Quaternion
	 * could represent a rotation, and the supplied vector a point to be rotated. 
	 * 
	 * @param x X component of vector to be transformed
	 * @param y Y component of vector to be transformed
	 * @param z Z component of vector to be transformed
	 * 
	 * This code is based on http://www.cprogramming.com/tutorial/3d/quaternions.html, accessed on 09/07/2013.
	 */	
	public Double3D transform(Double3D vec)
	{	return transform(vec.x, vec.y, vec.z);	}
	
	/** Applies the transformation represented by this Quaternion to the specified vector. For example, this Quaternion
	 * could represent a rotation, and the supplied vector a point to be rotated. 
	 * 
	 * @param x X component of vector to be transformed
	 * @param y Y component of vector to be transformed
	 * @param z Z component of vector to be transformed
	 * 
	 * This code is based on http://www.cprogramming.com/tutorial/3d/quaternions.html, accessed on 09/07/2013.
	 */	
	public Double3D transform(double vecX, double vecY, double vecZ)
	{
		double new_x =   w*w*vecX + 2*y*w*vecZ - 2*z*w*vecY +   x*x*vecX + 2*y*x*vecY + 2*z*x*vecZ -   z*z*vecX - y*y*vecX;
		double new_y = 2*x*y*vecX +   y*y*vecY + 2*z*y*vecZ + 2*w*z*vecX -   z*z*vecY +   w*w*vecY - 2*x*w*vecZ - x*x*vecY;
		double new_z = 2*x*z*vecX + 2*y*z*vecY +   z*z*vecZ - 2*w*y*vecX -   y*y*vecZ + 2*w*x*vecY -   x*x*vecZ + w*w*vecZ;
		return new Double3D(new_x, new_y, new_z);		
	}
	
	public Quaternion clone()
	{
		return new Quaternion(this.w, this.x, this.y, this.z);
	}
	
	
	/** q1 can be non-normalised quaternion
	 * 
	 *  @return 3-item array of doubles, containing the yaw, pitch and roll Euler angles.  
	 */
	public double[] toEulerAngles() 
	{
		double yaw = 0.;
		double pitch = 0.;
		double roll = 0.;
		
	    double sqw = w*w;
	    double sqx = x*x;
	    double sqy = y*y;
	    double sqz = z*z;
		double unit = sqx + sqy + sqz + sqw; // if normalised is one, otherwise is correction factor
		double test = x*y + z*w;		
		if (test > 0.499*unit)				// singularity at north pole 
		{ 	
			yaw = 2 * Math.atan2(x,w);
			roll = Math.PI/2;
			roll = 0;
		} else if (test < -0.499*unit) { 	// singularity at south pole
			yaw = -2 * Math.atan2(x,w);
			pitch = -Math.PI/2;
			roll = 0;
		} else {
			yaw = Math.atan2(2*y*w - 2*x*z , sqx - sqy - sqz + sqw);
			pitch = Math.asin(2*test/unit);
			roll = Math.atan2(2*x*w - 2*y*z , -sqx + sqy - sqz + sqw);
		}
		return new double[]{yaw,pitch,roll};
	}
	
	/** Converts this quaternion into axis & angle representation. Axis & angle defines a rotation as an angle (in 
	 *  radians) around an axis. Note that this is not the same as quaternion representation, thought quaternions are
	 *  also expressed in 4D space. 
	 * 
	 * @return 4-item array containing angle, and axis components of X Y and Z. 
	 */
	public double[] toAxisAngle()
	{
		double angle = 2. * Math.acos(w);
		double axisX = x / Math.sqrt(1. - w*w);
		double axisY = y / Math.sqrt(1. - w*w);
		double axisZ = z / Math.sqrt(1. - z*z);
		return new double[]{angle, axisX, axisY, axisZ};
	}
	
	
	public String toString()
	{
		return "Quaternion[" +w+ ", " +x+ "i, " +y+ "j, " +z+ "k]";
	}
	
	/* Used for debugging */
	public boolean fatalError_NAN()
	{	return Double.isNaN(x);		}
	
	/** Provided for testing the Quanternion class performance with. Also provides examples of how to use this class. 
	 */
	public static void main(String[] args)
	{
		Quaternion orientation = Quaternion.identity();			// ongoing record of orientation as changes are made.
		Double3D heading_start = new Double3D(0, 1, 0);			// when transformed through the above, this should point along current orientation.   
		
		Quaternion rotationQ = Quaternion.representRotation(Math.PI/2, 1, 0, 0);	// rotate 90 degrees around z axis.
		Double3D heading = new Double3D(heading_start.x, heading_start.y, heading_start.z);			// iteratively updated. Start off pointing along the x axis.
		
		// used for quaternion maths, rather than matrix.		 
		Quaternion heading_q = new Quaternion(0, heading_start.x, heading_start.y, heading_start.z);	 
		
		System.out.println("\nstart ");
		System.out.println("heading = " + heading);
		System.out.println("headingT = " + orientation.transform(heading_start));
		
		for (int i = 0; i < 4; i++)
		{
			heading = rotationQ.transform(heading); 			// move the heading around. 
			orientation = orientation.multiply(rotationQ);
			heading_q = rotationQ.multiply(heading_q).multiply(rotationQ.conjugate());
			
			System.out.println("\niteration "+i);
			System.out.println("heading = " + heading);
			System.out.println("headingT = " + orientation.transform(heading_start));			
			System.out.println("headingQ = " + heading_q);
		}		
	}
}
