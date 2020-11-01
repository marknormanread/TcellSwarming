package movement;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.MigratoryCell;
import core.Simulation;
import sim.util.Double3D;
import utils.Quaternion;

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
public class InverseHomogeneousCRW implements Translation, Orientation
{

	private static double speedMean;				// Microns per minute.
	private static double speedStD;					// Standard deviation.	
	// controls change of direction in random walk. Given in PI-radians (1 PI-radian = 180 degrees).
	static double pitchRateMean;
	static double pitchRateStd;
	static double rollRateMean;
	static double rollRateStd;
	
	// This relates the speed to the magnitude of pitch change. 
	protected static double correlationFactor;
	protected static double attackFactor;
	
	private double currentSpeed;
	
	public static InverseHomogeneousCRW instance = new InverseHomogeneousCRW();
	private InverseHomogeneousCRW() 
	{}
	
	private void calculateSpeed()
	{		
		/* apply movement to the cell in the direction that it faces */		
		this.currentSpeed = (Simulation.instance.random.nextGaussian() * speedStD) + speedMean;
	}
	
	
	public Double3D move(Quaternion orientation)
	{	
		// cell moves along it's x axis. Find its orientation in absolute space, by transforming x-axis. This gives a 
		// unit vector poining in the direction of the cell's orientation. 
		Double3D facing = orientation.transform(MigratoryCell.x_axis);		
 
		double dist = currentSpeed * Simulation.timeSlice_min;
		// translate would-be backwards movement into forwards. 
		dist = Math.abs(dist);
		// convert unit vector describing cell's orientation in absolute space to a move forward. 
		Double3D move = new Double3D(facing.x * dist, facing.y * dist, facing.z * dist);		
		return move;		
	}
	
	/** 
	 * Provides a new orientation of the cell, based on its current orientation (supplied as arg), and 
	 * distributions describing their propensity to chang direction.
	 * 
	 * Cells are rotated along their axis (direction in which they face, defined as the x-axis relative to the 
	 * cell), essentially rolling them, after which their pitch is changed (rotate around y-axis). 
	 * 
	 * @return The new orientation. 
	 */
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell)
	{	
		// must know how fast the cell is going to move prior to scaling pitch angle appropriately.
		calculateSpeed();		
		/* Note the order that quaternion multiplication is done here - the rotations are expressed relative to the
		 * cell, hence orientation is multiplied by the rotation quaternion.  */
		// roll the cell along it's x-axis (axis in which it faces). 
		// This rolls the cell, but doesn't change it's heading or pitch. 
		double roll;
		if (rollRateMean < 0.0) {
			// if mean roll rate is negative, assume this indicates a uniform distribution
			// should be used. 
			roll = Simulation.instance.random.nextDouble() * 2.0 * Math.PI;
		} else {
			roll = (Simulation.instance.random.nextGaussian() * rollRateStd) + rollRateMean;
			// randomly invert roll direction. Avoids corkscrewing. 
			if (Simulation.instance.random.nextBoolean())
				roll *= -1.0;		// cells can roll in either direction.
			roll *= Simulation.timeSlice_min;
		}	
		// roll as a quaternion.
		Quaternion rotateQ = Quaternion.representRotation
				(roll, MigratoryCell.x_axis.x, MigratoryCell.x_axis.y, MigratoryCell.x_axis.z); 
		// multiply orientation by rotateQ, because rotateQ is calculated relative to cell, not in absolute space. 
		orientation = orientation.multiply(rotateQ).normalise();	// alter the cell's orientation. 
		
		// change cell pitch (roll along the y axis). Pitch can be changed in both positive and negative directions.
		double pitch = (Simulation.instance.random.nextGaussian() * pitchRateStd) + pitchRateMean;
		
		// Math.pow (used below) can't take negative base numbers. If speed < 0, then take its absolute 
		double power_calc_speed = this.currentSpeed >=0 ? this.currentSpeed : Math.abs(this.currentSpeed); 
		/* This is the critical bit, relating the last speed to the pitch. 
		 * Correlation factor describes the sensitivity to speed. When correlationFactor=0, speed has no influence
		 * on pitch. When correlationFactor=1, speed as a proportion of maximum speed completely scales pitch.  
		 */
		double maxSpeed = 30.0;
		double scaledMaxSpeed = Math.pow(maxSpeed, attackFactor);
		double scaledSpeed = Math.pow(power_calc_speed, attackFactor);  // avoid supplying negative speeds here
		pitch *= (1.0 - correlationFactor) + (correlationFactor * (scaledMaxSpeed - scaledSpeed) / scaledMaxSpeed);
		pitch *= Simulation.timeSlice_min;		// account for timestep.		
		Quaternion pitchQ = Quaternion.representRotation
				(pitch, MigratoryCell.y_axis.x, MigratoryCell.y_axis.y, MigratoryCell.y_axis.z);
		// multiply orientation by rotateQ, because pitchQ is calculated relative to cell, not in absolute space.
		orientation = orientation.multiply(pitchQ).normalise();		
		return orientation;
	}
	
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;
		n = (Node) xPath.compile("/params/Motility/InverseHomogeneousCRW/speedMean")
				.evaluate(params, XPathConstants.NODE);
		speedMean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHomogeneousCRW/speedStD")
				.evaluate(params, XPathConstants.NODE);
		speedStD = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Motility/InverseHomogeneousCRW/pitchRateMean")
				.evaluate(params, XPathConstants.NODE);
		pitchRateMean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHomogeneousCRW/pitchRateStd")
				.evaluate(params, XPathConstants.NODE);
		pitchRateStd = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHomogeneousCRW/rollRateMean")
				.evaluate(params, XPathConstants.NODE);
		rollRateMean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHomogeneousCRW/rollRateStd")
				.evaluate(params, XPathConstants.NODE);
		rollRateStd = Double.parseDouble(n.getTextContent());	
		
		n = (Node) xPath.compile("/params/Motility/InverseHomogeneousCRW/corr")
				.evaluate(params, XPathConstants.NODE);
		correlationFactor = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHomogeneousCRW/attack")
				.evaluate(params, XPathConstants.NODE);
		attackFactor = Double.parseDouble(n.getTextContent());
	}
}
