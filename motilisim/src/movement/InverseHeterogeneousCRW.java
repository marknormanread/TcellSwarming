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
import utils.Utils;

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
public class InverseHeterogeneousCRW implements Orientation, Translation
{
	// describe a distribution from which the mean of a track's speed distribution is selected. 
	protected static double speedM_Mean = 0.0;
	protected static double speedM_StD = 0.0;
	// describe a distribution from which the standard deviation of a track's speed distribution is selected. 
	protected static double speedS_Mean = 0.0;
	protected static double speedS_StD = 0.0;
		
	protected static double pitchM_Mean;
	protected static double pitchM_StD;
	protected static double pitchS_Mean;
	protected static double pitchS_StD;

	protected static double rollM_Mean;
	protected static double rollM_StD;
	protected static double rollS_Mean;
	protected static double rollS_StD;
	
	
	/* Create custom random walk mechanism for a cell. Mean and std of this cell's movememnt dynamics are
	 * themselves drawn from a distribution. */		
	protected double speedMean = 0.0;	// instance distribution parameters. Accessible to subclasses. 
	protected double speedStD = 0.0;
	
	protected double pitchRateMean;
	protected double pitchRateStD;
	protected double rollRateMean;
	protected double rollRateStD;
		
	// This relates the speed to the magnitude of pitch change. 
	protected static double correlationFactor;
	protected static double attackFactor;
	
	private double currentSpeed;
	
	public InverseHeterogeneousCRW()
	{
		speedMean = (Simulation.instance.random.nextGaussian() * speedM_StD) + speedM_Mean;
		// invert negative values. 
		if (speedMean < 0.0) 	speedMean *= -1.0;
		speedStD = (Simulation.instance.random.nextGaussian() * speedS_StD) + speedS_Mean;
		if (speedStD < 0.0) speedStD *= -1.0;		// invert negative values		
		
		pitchRateMean = (Simulation.instance.random.nextGaussian() * pitchM_StD) + pitchM_Mean;
		pitchRateMean = Math.abs(pitchRateMean);	// invert negative values.
		pitchRateStD = (Simulation.instance.random.nextGaussian() * pitchS_StD) + pitchS_Mean;
		pitchRateStD = Math.abs(pitchRateStD);
		
		rollRateMean = (Simulation.instance.random.nextGaussian() * rollM_StD) + rollM_Mean;
		rollRateMean = Math.abs(rollRateMean);	// invert negative values.
		rollRateStD = (Simulation.instance.random.nextGaussian() * rollS_StD) + rollS_Mean;
		rollRateStD = Math.abs(rollRateStD);			
	}

	
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
	
	@Override
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell) 
	{
		// must know how fast the cell is going to move prior to scaling pitch angle appropriately.
		calculateSpeed();		
		/* Note the order that quaternion multiplication is done here - the rotations are expressed relative to the
		 * cell, hence orientation is multiplied by the rotation quaternion.  */
		// roll the cell along it's x-axis (axis in which it faces). 
		// This rolls the cell, but doesn't change it's heading or pitch. 
		double roll;
		if (rollM_Mean < 0.0) {
			// if mean roll rate has been set to a negative value, assume this indicates a uniform distribution
			// should be used. 
			roll = Simulation.instance.random.nextDouble() * 2.0 * Math.PI;
		} else {
			roll = (Simulation.instance.random.nextGaussian() * rollRateStD) + rollRateMean;
			// randomly invert roll direction to prevent corkscrewing. 
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
		double pitch = (Simulation.instance.random.nextGaussian() * pitchRateStD) + pitchRateMean;

		// Math.pow (used below) can't take negative base numbers. If speed < 0, then take its absolute 
		double power_calc_speed = this.currentSpeed >=0 ? this.currentSpeed : Math.abs(this.currentSpeed); 
		/* This is the critical bit, relating the last speed to the pitch. 
		 * Correlation factor describes the sensitivity to speed. When correlationFactor=0, speed has no influence
		 * on pitch. When correlationFactor=1, speed as a proportion of maximum speed completely scales pitch.  
		 */
		double maxSpeed = 30.0;
		double scaledMaxSpeed = Math.pow(maxSpeed, attackFactor);
		double scaledSpeed = Math.pow(power_calc_speed, attackFactor);
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
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/speedM_Mean")
				.evaluate(params, XPathConstants.NODE);
		speedM_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/speedM_StD")
				.evaluate(params, XPathConstants.NODE);
		speedM_StD = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/speedS_Mean")
				.evaluate(params, XPathConstants.NODE);
		speedS_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/speedS_StD")
				.evaluate(params, XPathConstants.NODE);
		speedS_StD = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/pitchM_Mean")
				.evaluate(params, XPathConstants.NODE);
		pitchM_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/pitchM_StD")
				.evaluate(params, XPathConstants.NODE);
		pitchM_StD = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/pitchS_Mean")
				.evaluate(params, XPathConstants.NODE);
		pitchS_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/pitchS_StD")
				.evaluate(params, XPathConstants.NODE);
		pitchS_StD = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/rollM_Mean")
				.evaluate(params, XPathConstants.NODE);
		rollM_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/rollM_StD")
				.evaluate(params, XPathConstants.NODE);
		rollM_StD = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/rollS_Mean")
				.evaluate(params, XPathConstants.NODE);
		rollS_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/rollS_StD")
				.evaluate(params, XPathConstants.NODE);
		rollS_StD = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/corr")
				.evaluate(params, XPathConstants.NODE);
		correlationFactor = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/InverseHeterogeneousCRW/attack")
				.evaluate(params, XPathConstants.NODE);
		attackFactor = Double.parseDouble(n.getTextContent());
	}
}
