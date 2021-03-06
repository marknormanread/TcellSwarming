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
    
 * Cells follow a ballistic trajectory. This is used primarily for testing purposes. 
 * 
 * @author Mark N. Read
 *
 */
public class Ballistic implements Orientation, Translation
{
	protected static double speed;
	
	// singleton pattern. 
	public static Ballistic instance = new Ballistic();
	private Ballistic()
	{}
	
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell)
	{
		return orientation;	
	}
	
	@Override
	public Double3D move(Quaternion orientation)
	{	
		// cell moves along it's x axis. Find its orientation in absolute space, by transforming x-axis. This gives a 
		// unit vector pointing in the direction of the cell's orientation. 
		Double3D facing = orientation.transform(MigratoryCell.x_axis);		
		/* apply movement to the cell in the direction that it faces.  */				 
		double dist = speed * Simulation.timeSlice_min; // units in um/min.
		// translate would-be backwards movement into forwards. 
		dist = Math.abs(dist);
		// convert unit vector describing cell's orientation in absolute space to a move forward. 
		Double3D move = new Double3D(facing.x * dist, facing.y * dist, facing.z * dist);			
		return move;		
	}
	
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;
		
		n = (Node) xPath.compile("/params/Motility/Ballistic/speed").evaluate(params, XPathConstants.NODE);
		speed = Double.parseDouble(n.getTextContent());
	}
}
