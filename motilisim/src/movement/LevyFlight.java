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
import utils.LevyDistribution;
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
 * 
 * Implements a Levy Flight random walk. Cells pick a random orientation, and then move in that direction by some
 * distance, determined from a long-tailed distribution where the long tail points to positive infinity.   
 */
public class LevyFlight implements Orientation, Translation
{	
	
	/* Note that the Mu values here must reside within the range 1..3 (inclusive), or the Levy
	 * calculation returns NaN. 
	 */
	private static double restMu = Double.NaN;
	private static double restScale = Double.NaN;
	
	private static double motileMu = Double.NaN; 
	private static double motileScale = Double.NaN;
	
	private static double speedMu = Double.NaN;
	private static double speedScale = Double.NaN;
	
	private enum State
	{
		REST,
		MOTILE
	}
	private State state = State.REST;

	
	private Quaternion newOrientation;
	private double currentSpeed;
	
	private double endMotileTime = -1.0;	// absolute simulation time at which cell will start moving again. 
	private double endRestTime = -1.0;    // absolute simulation time at which cell will pause.
	
	public LevyFlight()
	{
		plan();
	}

	private void plan()
	{	
		// execute this twice so that cases where no rest (or motility, though that would be daft) are handled.
		// Otherwise there will be a time step where nothing happens, even if the rest/motile duration is zero.
		for (int i = 0; i < 2; i++)
		{
			if (this.state == State.REST)
			{
				if (Simulation.instance.schedule.getTime() >= endRestTime)
				{
					// timer has elapsed, start moving again.
					newOrientation = Quaternion.randomUniform();
					currentSpeed = LevyDistribution.sample_positive(speedMu, speedScale);
					// Levy distributions can supply values far, far larger than plausible as a cell's speed. 
					// The distances entailed by these speeds can cause computational headaches for looking up 
					// potential collisions. Hence, we limit the speed to the extreme upper end of what we would ever
					// expect to see (and probably only cells in the vasculature at that). 
					currentSpeed = Math.min(currentSpeed, 100.);
					
					double duration = LevyDistribution.sample_positive(motileMu, motileScale);
					endMotileTime = Simulation.instance.schedule.getTime() + duration;
					this.state = State.MOTILE;
				}
			} 
			if (this.state == State.MOTILE) 
			{
				// cell is moving.
				if (Simulation.instance.schedule.getTime() >= endMotileTime)
				{	
					// time to stop moving. 				
					currentSpeed = 0.0;
					double duration = 0.0;
					if (restMu > 0.0)	// if rest time is being used. 
						duration = LevyDistribution.sample_positive(restMu, restScale);					
					endRestTime = Simulation.instance.schedule.getTime() + duration;
					this.state=State.REST;
				}
			}
		}
	}
	
	@Override
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell) 
	{
		plan();
		return newOrientation;		
	}

	@Override
	public Double3D move(Quaternion orientation) 
	{
		// Levy flight draws step size from a long tailed distribution. In this case log-normal.		
		double length = currentSpeed * Simulation.timeSlice_min;
		// cell moves along it's x axis. Find its orientation in absolute space, by transforming x-axis. This gives a 
		// unit vector pointing in the direction of the cell's orientation.
		Double3D facing = orientation.transform(MigratoryCell.x_axis);
		// convert unit vector describing cell's orientation in absolute space to a move forward.
		Double3D move = new Double3D(facing.x*length, facing.y*length, facing.z*length);
		return move;
	}

	
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;		
		
		n = (Node) xPath.compile("/params/Motility/LevyFlight/restMu").evaluate(params, XPathConstants.NODE);
		restMu = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/LevyFlight/restScale").evaluate(params, XPathConstants.NODE);
		restScale = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Motility/LevyFlight/motileMu").evaluate(params, XPathConstants.NODE);
		motileMu = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/LevyFlight/motileScale").evaluate(params, XPathConstants.NODE);
		motileScale = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Motility/LevyFlight/speedMu").evaluate(params, XPathConstants.NODE);
		speedMu = Double.parseDouble(n.getTextContent());	
		n = (Node) xPath.compile("/params/Motility/LevyFlight/speedScale").evaluate(params, XPathConstants.NODE);
		speedScale = Double.parseDouble(n.getTextContent());	
	}
}	

