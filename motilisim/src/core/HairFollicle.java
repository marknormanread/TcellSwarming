package core;

import sim.engine.SimState;
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
 * Class implements cell so that the Compartment3D class can access diameter and radius information. 
 * @author Mark N. Read
 *
 */
public class HairFollicle extends Cell
{
	private double diameter;
	
	public HairFollicle()
	{		
		// diameters are between 40 and 60 micrometers. 
		diameter = Simulation.instance.random.nextDouble() * 20 + 40;
	}
	
	public HairFollicle(double diameter, Double3D location)
	{
		this.diameter = diameter; 
		this.location = location;
	}
	
	public void setLocation(Double3D location)
	{	this.location = location;	}
	
	public double getDiameter()
	{	return diameter;	}
	
	public double getRadius()
	{	return diameter / 2.0;	}
	
	public Double3D getCentreLocation()
	{	return location;	}

	@Override
	public void step(SimState state) 
	{ /* do nothing, hair follicles do not change state */ }
	
	public String getType()
	{	return "folicle";	}
}
