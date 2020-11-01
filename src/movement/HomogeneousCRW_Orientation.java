package movement;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.MigratoryCell;
import core.Simulation;
import utils.Quaternion;

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
 * Changes in cell orientation achieved changes in roll and pitch. The magnitude of these changes are drawn from
 * normal distributions that are common to all cells.  
 */
public class HomogeneousCRW_Orientation implements Orientation
{
	// controls change of direction in random walk. Given in PI-radians (1 PI-radian = 180 degrees).
	public static double pitchRateMean;
	public static double pitchRateStd;
	public static double rollRateMean;
	public static double rollRateStd;
	
	/* Singleton pattern */
	public static HomogeneousCRW_Orientation instance = new HomogeneousCRW_Orientation();
	protected HomogeneousCRW_Orientation()  // visible only to subclasses
	{}
	
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
		n = (Node) xPath.compile("/params/Motility/HomogeneousCRW_Orientation/pitchRateMean")
				.evaluate(params, XPathConstants.NODE);
		pitchRateMean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HomogeneousCRW_Orientation/pitchRateStd")
				.evaluate(params, XPathConstants.NODE);
		pitchRateStd = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HomogeneousCRW_Orientation/rollRateMean")
				.evaluate(params, XPathConstants.NODE);
		rollRateMean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HomogeneousCRW_Orientation/rollRateStd")
				.evaluate(params, XPathConstants.NODE);
		rollRateStd = Double.parseDouble(n.getTextContent());		
	}
}