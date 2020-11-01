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
public class HeterogeneousCRW_Orientation implements Orientation
{
	protected double pitchRateMean;
	protected double pitchRateStD;
	protected double rollRateMean;
	protected double rollRateStD;
	
	protected static double pitchM_Mean;
	protected static double pitchM_StD;
	protected static double pitchS_Mean;
	protected static double pitchS_StD;

	protected static double rollM_Mean;
	protected static double rollM_StD;
	protected static double rollS_Mean;
	protected static double rollS_StD;
	
	public HeterogeneousCRW_Orientation()
	{
		pitchRateMean = (Simulation.instance.random.nextGaussian() * pitchM_StD) + pitchM_Mean;
		pitchRateMean = Math.abs(pitchRateMean);	// invert negative values.
		pitchRateStD = (Simulation.instance.random.nextGaussian() * pitchS_StD) + pitchS_Mean;
		pitchRateStD = Math.abs(pitchRateStD);
		
		rollRateMean = (Simulation.instance.random.nextGaussian() * rollM_StD) + rollM_Mean;
		rollRateMean = Math.abs(rollRateMean);	// invert negative values.
		rollRateStD = (Simulation.instance.random.nextGaussian() * rollS_StD) + rollS_Mean;
		rollRateStD = Math.abs(rollRateStD);	
	}

	@Override
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell) 
	{
		/* Note the order that quaternion multiplication is done here - the rotations are expressed relative to the
		 * cell, hence orientation is multiplied by the rotation quaternion.  */
		// roll the cell along it's x-axis (axis in which it faces). 
		// This rolls the cell, but doesn't change it's heading or pitch. 
		double roll;
		if (rollM_Mean < 0.0)
			// if mean roll rate has been set to a negative value, assume this indicates a uniform distribution
			// should be used. 
			roll = Simulation.instance.random.nextDouble() * 2.0 * Math.PI;
		else{
			roll = (Simulation.instance.random.nextGaussian() * rollRateStD) + rollRateMean;
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
		// randomly invert values. This makes no difference to zero-mean distributions, and avoids a bias in 
		// non-zero-mean distributions. Hence, applicable to both. 
		if (Simulation.instance.random.nextBoolean())
			pitch *= -1.0;
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
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Orientation/pitchM_Mean")
				.evaluate(params, XPathConstants.NODE);
		pitchM_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Orientation/pitchM_StD")
				.evaluate(params, XPathConstants.NODE);
		pitchM_StD = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Orientation/pitchS_Mean")
				.evaluate(params, XPathConstants.NODE);
		pitchS_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Orientation/pitchS_StD")
				.evaluate(params, XPathConstants.NODE);
		pitchS_StD = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Orientation/rollM_Mean")
				.evaluate(params, XPathConstants.NODE);
		rollM_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Orientation/rollM_StD")
				.evaluate(params, XPathConstants.NODE);
		rollM_StD = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Orientation/rollS_Mean")
				.evaluate(params, XPathConstants.NODE);
		rollS_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Orientation/rollS_StD")
				.evaluate(params, XPathConstants.NODE);
		rollS_StD = Double.parseDouble(n.getTextContent());
	}
}
