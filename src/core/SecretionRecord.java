package core;

import sim.util.Double3D;

/**
 * Groups together information pertaining to where an entity secreted a particular
 * quantity of soluble factor, at a particular time, at a particular location. 
 * @author Mark N. Read, 2017
 *
 */
public class SecretionRecord 
{
	private Double3D location;
	private double time;
	private double quantity;
	
	public SecretionRecord(Double3D loc, double time, double quantity)
	{
		this.location = loc;
		this.time = time;
		this.quantity = quantity;
	}
	
	public Double3D getLocation()
	{	return this.location;	}
	
	public double getTime()
	{	return this.time;	}
	
	public double getQuantity()
	{ 	return this.quantity;	}
}
