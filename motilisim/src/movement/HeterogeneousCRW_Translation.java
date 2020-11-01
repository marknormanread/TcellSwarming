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
public class HeterogeneousCRW_Translation implements Translation
{
	// describe a distribution from which the mean of a track's speed distribution is selected. 
	protected static double speedM_Mean = 0.0;
	protected static double speedM_StD = 0.0;
	// describe a distribution from which the standard deviation of a track's speed distribution is selected. 
	protected static double speedS_Mean = 0.0;
	protected static double speedS_StD = 0.0;
	
	/* Create custom random walk mechanism for a cell. Mean and std of this cell's movememnt dynamics are
	 * themselves drawn from a distribution. */		
	protected double speedMean = 0.0;	// instance distribution parameters. Accessible to subclasses. 
	protected double speedStD = 0.0;
	
	// accessible to others in this package. This was the last speed drawn from the cell's distribution.
	double currentSpeed;
	
	public HeterogeneousCRW_Translation()
	{
		speedMean = (Simulation.instance.random.nextGaussian() * speedM_StD) + speedM_Mean;
		// invert negative values. 
		if (speedMean < 0.0) 	speedMean *= -1.0;
		speedStD = (Simulation.instance.random.nextGaussian() * speedS_StD) + speedS_Mean;
		if (speedStD < 0.0) speedStD *= -1.0;		// invert negative values
	}
	
	public Double3D move(Quaternion orientation)
	{			
		// cell moves along it's x axis. Find its orientation in absolute space, by transforming x-axis. This gives a 
		// unit vector poining in the direction of the cell's orientation. 
		Double3D facing = orientation.transform(MigratoryCell.x_axis);		
		/* apply movement to the cell in the direction that it faces */		
		currentSpeed = (Simulation.instance.random.nextGaussian() * speedStD) + speedMean;
		double dist = currentSpeed * Simulation.timeSlice_min;
		// translate would-be backwards movement into forwards. 
		// convert unit vector describing cell's orientation in absolute space to a move forward. 
		Double3D move = new Double3D(facing.x * dist, facing.y * dist, facing.z * dist);			
		return move;	
	}
	
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Translation/speedM_Mean")
				.evaluate(params, XPathConstants.NODE);
		speedM_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Translation/speedM_StD")
				.evaluate(params, XPathConstants.NODE);
		speedM_StD = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Translation/speedS_Mean")
				.evaluate(params, XPathConstants.NODE);
		speedS_Mean = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/HeterogeneousCRW_Translation/speedS_StD")
				.evaluate(params, XPathConstants.NODE);
		speedS_StD = Double.parseDouble(n.getTextContent());
	}
}
