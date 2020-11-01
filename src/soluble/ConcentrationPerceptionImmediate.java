package soluble;

import java.util.Arrays;
import java.util.Collections;

import core.Simulation;
import sim.util.Double3D;

/** 
 * Encapsulates concentrations measured around a cell. Intended to represent 6 points around a cell.
 * 
 * @author Mark N. Read, 2018. 
 */
public class ConcentrationPerceptionImmediate implements ConcentrationPerception
{
	final double sampleTime;
	
	private double front;
	private double back;
	private double right;
	private double left;
	private double up;
	private double down;
	
	private double sum; 
	private double max;
	
	private Double3D gradient = null;
	private double maxGradientDifferential;
	
		
	/** Concentrations are calculated externally and supplied. */
	public ConcentrationPerceptionImmediate(
			double front, double back, double right, double left, double up, double down, double sampleTime)
	{
		this.sampleTime = sampleTime;
		this.front = front;
		this.back = back;
		this.right = right;
		this.left = left;
		this.up = up;
		this.down = down;
				
		sum = front + back + right + left + up + down;
		max = calculateMax();
		calculateGradient();
	}
	
	/** Find maximum absolute concentration in this perception */
	public double maxConcentration()
	{
		return this.max;
	}
	
	/** Find maximum absolute concentration in this perception */
	private double calculateMax()
	{
		double max = front;
		final double[] tmp = {back, left, right, up, down};
		for (double t : tmp)
			max = Math.max(max, t);	
		return max;
	}
	
	/** Returns a Double3D representing the direction of the gradient relative to the orientation of the concentrations
	 *  supplied. Ie, if the concentrations represent points around the cell, then the vector will be relative to the
	 *  cell's orientation. 
	 *  
	 *  The vector is not a unit vector.  
	 *  It may have magnitude zero. 
	 *  Its magnitude corresponds to the gradient differential as measured around the cell.
	 */
	public Double3D gradientDirection()
	{
		if (this.gradient == null)
			calculateGradient();
		
		return this.gradient;  // Direction of gradient relative to the cell.		
	}	
	
	public double maxConcentrationDifferential()
	{
		if (this.gradient == null)
			calculateGradient();
		
		return this.maxGradientDifferential;
	}
	
	private void calculateGradient()
	{
		// Differential of gradients across orthogonal axes (relative to the cell, cell faces 'x'). 
		final double dx = this.front - this.back;	
		final double dy = this.left - this.right; 	
		final double dz = this.up - this.down;
		this.gradient = new Double3D(dx, dy, dz);
		this.maxGradientDifferential = Math.max(Math.max(dx, dy), dz);
	}
	
	/** If the time at which this peception was taken is in the past, the it is out of date; expired */
	public boolean expired()
	{	return Simulation.instance.schedule.getTime() > this.sampleTime;	}
}
