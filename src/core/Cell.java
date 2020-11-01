package core;

import sim.engine.Schedule;
import sim.engine.Steppable;
import sim.engine.Stoppable;
import sim.util.Double3D;
import utils.Quaternion;

/**
 *  This program is free software: you can redistribute it and/or modify
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
 * @author Mark N. Read
 *
 */
public abstract class Cell implements Steppable
{
	// Cells maintain an up-to-date copy of their location. 
	// It is expensive to keep asking the spatial compartment where a cell is, and this is done frequently. 
	protected Double3D location;
	
	private int id;  // Corresponds to track ID. 
	
	private static int id_counter = 0;
	
	protected Stoppable unSchedule;
	
	/* `orientation` represents the cell's current orientation in space, compared to where it started from. This 
	 * variable is updated every time the cell changes its orientation. Hence, it provides a conversion between space
	 * relative to the cell, and absolute space. 
	 * 
	 * orientation.transform(x_axis) will give a unit vector in absolute space that points along the x-axis relative to 
	 * the cell. Cells always move along their x_axis, so this is used to update a cells movement. 
	 * 
	 * The inverse can be used to convert coordinates relative to the cell back into absolute space. */
	protected Quaternion orientation = Quaternion.randomUniform();	 	
	public static final Double3D x_axis = new Double3D(1, 0, 0);			
	public static final Double3D y_axis = new Double3D(0, 1, 0);
	public static final Double3D z_axis = new Double3D(0, 0, 1);	
	
	/* Standard rotations to locate points around the cell that are frequently used, represented as quaternions. 
	 * No need to destroy and recreate objects, as the values don't change. Hence, static final objects. 
	 * Remember right hand skew. 
	 * Hence negative and positive values for z-axis of rotation. */
	final static Quaternion rotFront 	= Quaternion.representRotation(0, 1, 0, 0);
	final static Quaternion rotBack 	= Quaternion.representRotation(Math.PI,   0,  0, -1);
	final static Quaternion rotRight 	= Quaternion.representRotation(Math.PI/2, 0,  0, -1);
	final static Quaternion rotLeft 	= Quaternion.representRotation(Math.PI/2, 0,  0,  1);
	final static Quaternion rotUp 		= Quaternion.representRotation(Math.PI/2, 0, -1,  0);
	final static Quaternion rotDown 	= Quaternion.representRotation(Math.PI/2, 0,  1,  0);
	
	public Cell(){
		this.id = id_counter;
		id_counter++;
	}
	
	public Cell(Schedule sched)
	{
		this.id = id_counter;
		id_counter++;
		// Start time at which the new cell will be scheduled.
		double startTime = sched.getTime();
		if (sched.getTime() < Schedule.EPOCH)	// Time is -1.0 when simulation has not been started yet. 
			startTime = Schedule.EPOCH;		
		unSchedule = sched.scheduleRepeating(startTime, Simulation.cellOrdering, this, Simulation.timeSlice_min);
	}

	/** Should the cell be visualized and have it's location recorded? Defaults to true, but can be overridden. */
	public boolean recordable = true;
	
	public abstract double getDiameter();
	
	public abstract double getRadius();
	
	public Double3D getCurrentLocation()
	{	
		if(location != null)
			return location;
		else
			return Simulation.space.getCellLocation(this);
	}
	
	public void setCurrentLocation(Double3D loc)
	{	
		this.location = loc;
	}
	
	public int getID()
	{
		return this.id;
	}
	
	public abstract String getType();
	
	/** Retrieves Double3D objects representing the front, back, right, left, top and bottom sides of the cell, oriented
	 *  relative to the cell. Hence, the locations returned are points in absolute space (not relative to cell). 
	 */
	protected SamplePoints getSamplePoints()
	{
		SamplePoints sp = new SamplePoints();				
		final Double3D surf = x_axis.multiply(getRadius());	// Point directly ahead, on cell surface. 
		/* Cells have arbitrary orientations in space, they have a notion of front, back, left, right, up and down.
		 * Retrieve the points in absolute space that correspond to these locations relative to the cell. 
		 * The conversion between cell-space and absolute space assumes cell's face along the x axis. 
		 * Hence, first step is to take the rotations that express e.g. "up" relative to the cell (facing x axis), 
		 * and combine these with the cell's "orientation" quaternion, which expresses the cell's rotation w.r.t. 
		 * absolute space.
		 * Transform by surface, which will give an actual point on the cell's surface, in absolute space. 
		 * E.g., if a cell has radius 10, and we wish to find its "back" point in aboslute space, first we have to 
		 * take the rotation that says "'back' is 180 degrees behind the direction you face (x axis). We combine that
		 * with the cell's "orientation" to account for the face that the cell is not in fact facing the xaxis in 
		 * absolute space. Thereafter, we do a translational transformation by the cell's radius to find the location.  
		 * Lastly we add the cell's "location", to transform these points to absolute space at the cell's location.  
		 */
		final Double3D absOrientFront 	= orientation.multiply(rotFront).transform(surf);
		final Double3D absOrientBack 	= orientation.multiply(rotBack).transform(surf);
		final Double3D absOrientRight	= orientation.multiply(rotRight).transform(surf);
		final Double3D absOrientLeft 	= orientation.multiply(rotLeft).transform(surf);
		final Double3D absOrientUp 		= orientation.multiply(rotUp).transform(surf);
		final Double3D absOrientDown 	= orientation.multiply(rotDown).transform(surf);
		
		sp.front 	= absOrientFront.add(location);
		sp.back 	= absOrientBack.add(location);
		sp.right 	= absOrientRight.add(location);
		sp.left 	= absOrientLeft.add(location);
		sp.up 		= absOrientUp.add(location);
		sp.down 	= absOrientDown.add(location);
		return sp;
	}
	
	/** Convenient encapsulator of location objects around a cell. Intended to represent 6 points around a cell. 
	 * 
	 * @author Mark Read
	 */
	public class SamplePoints
	{
		public Double3D front;
		public Double3D back;
		public Double3D right;
		public Double3D left;
		public Double3D up;
		public Double3D down;		
	}
}
