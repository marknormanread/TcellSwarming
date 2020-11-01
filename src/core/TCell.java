package core;

import java.util.Vector;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Node;

import environment.BoundedCylinder;
import environment.Compartment3D.MoveResults;
import loggers.CellLogger;
import movement.ChemotaxisResponsive;
import sim.engine.Schedule;
import sim.engine.SimState;
import sim.util.Double3D;
import soluble.AttractantSecretor;
import soluble.ConcentrationPerception;
import soluble.ConcentrationPerceptionDeferred;
import soluble.ConcentrationPerceptionGradientMaxSupplied;
import soluble.ConcentrationPerceptionImmediate;

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

 * @author Mark N. Read, 2017
 *
 */

public class TCell extends MigratoryCell implements AttractantSecretor, ChemotaxisResponsive
{
	public static double diameter = 12.0;  // In microns.

	private static double attractantSecretionRate = 0.0;  // Per minute.
	protected static double relayActivationAttractionThreshold = 100.;

	public static double chemotaxisConcentrationLocation = 0;
	public static double chemotaxisConcentrationShape = 0;
	public static double chemotaxisConcentrationDesensitizationLocation = 0;
	public static double chemotaxisConcentrationDesensitizationShape = 0;
	private static LogNormalDistribution chemokineResponseDistro;
	private static LogNormalDistribution chemokineDesensitizationDistro;
	// Cell-specific sensitivities. 
	// This is the chemokine concentration required to initialise chemotaxis.
	private double chemokineConcentrationResponseThreshold = Double.POSITIVE_INFINITY;  
	// This is the chemokine concentration at which chemotaxis ceases. 
	private double chemokineConcentrationDesensitizationThreshold = Double.POSITIVE_INFINITY;
	
		
	// Where and when the cell has been secreting
	private Vector<SecretionRecord> secretionHistory = new Vector<SecretionRecord>();
	private double secretionMemory = Double.POSITIVE_INFINITY;  // Period of recent secretion events maintained. 
	
	protected boolean secretingAttractant = false;

	// How far off the x-axis a cell can perceive a gradient. In PI-radians (1 PI-radian = 180 degrees).
	// Needed to stop cells 'chasing their tails' if also secreting themselves.
	protected static double perceptionCone = 0.;

	protected static enum SecretionModeTypes {
		None,          // Cells never secrete any factors.
		WithinBolus,   // Cells secrete only when within the bolus.
		// Not currently implemented. Need threshold for activation; cells can secrete outside of bolus, given 
		// sufficient exposure to chemokine. 
		Relay,         
		BolusCentre,   // No secretion, but all cells sense a gradient towards bolus centre. Testing. 
		// Secretion in bolus centre, PDE-based concentration detection, but gradient is directed to bolus centre.
		// PDE-based concentration detection is computationally costly. 
		// Calculating gradient requires sampling 6 points around cell. 
		// This samples only once, at the cell's front. 
		WithinBolusGradientCentre
	}
	protected static SecretionModeTypes secretionMode = SecretionModeTypes.None;

	private ConcentrationPerception attractantPerception;

	private CellLogger.Track logger;

	public boolean motile = true;  // Arrest cells in bolus

	protected TCell()
	{}

	public TCell(Schedule sched)
	{
		super(sched);
		if (Simulation.trackCells)
			logger = new CellLogger.Track(this);

		if (chemokineResponseDistro != null)
		{
			chemokineConcentrationResponseThreshold = chemokineResponseDistro.sample();
			chemokineConcentrationDesensitizationThreshold = chemokineDesensitizationDistro.sample();
		} else {
			chemokineConcentrationResponseThreshold = Double.POSITIVE_INFINITY;
		}
	}

	@Override
	public void step(SimState state)
	{
		location = Simulation.space.getCellLocation(this);
		
		if (this.motile)
		{
			// Special processing required if cell has hit something (non-zero bounce length)
			if (bounce.lengthSq() != 0.0)
			{
				// Change cell's orientation in response to hitting an obstacle in previous step
				bounce();
			} else {
				// Normal dynamics, nothing hit
				orientation = orientationActuator.newOrientation(orientation, this);

				// Handles case for NON bootstrapped chemotaxis. Ie, other forms of motility. 
				if (! (   secretionMode == SecretionModeTypes.None
					   || orientationParadigm == OrientationParadigm.BOOTSTRAP
					   || orientationParadigm == OrientationParadigm.BOOTSTRAP_GRADIENT
					   || orientationParadigm == OrientationParadigm.BOOTSTRAP_DYNAMIC))
				{
					rotateToVectorLimited(perceiveAttractant().gradientDirection(), perceptionCone);
				}
			}
			Double3D move = null;
			// What the cell would *like* to do. Environment may override or curtail. 
			move = translationActuator.move(orientation);  
			
			// Performs both inter-cellular and boundary collision detection.
			MoveResults mr = Simulation.space.moveCellCollisionDetection(this, move);
			location = mr.newLocation;
			bounce = mr.bounce;	 // Bounce off other cells that may have been contacted.
			
			if (Double.isNaN(location.x))  // This should never happen. Here to catch possible bugs. 
			{
				throw new RuntimeException("Set agent location to NaN : " + location);
			}
		}

		// If in relay mode, and if sufficient concentration of attractant, start secreting.
		if (secretionMode == SecretionModeTypes.Relay && perceiveAttractant().maxConcentration()
														 >= relayActivationAttractionThreshold)
			secretingAttractant = true;

		// If in bolus start secreting.
		if (Simulation.space.insideBolus(location)  
				&& (   secretionMode == SecretionModeTypes.WithinBolus 
				    || secretionMode == SecretionModeTypes.WithinBolusGradientCentre )
			)
			secretingAttractant = true;

		if (secretingAttractant)
			secrete((Simulation)state);
	}

	public void secrete(Simulation sim)
	{
		double secretionRate = attractantSecretionRate * Simulation.timeSlice_min;
		
		// Get location of back of cell, this is where secretion takes place.
		Double3D back = orientation.transform(-getRadius(), 0.0, 0.0).add(location);
		secretionHistory.add(new SecretionRecord(
				back,  // Location where secretion takes place. Back of cell.
				sim.schedule.getTime(),  // Time of secretion
				secretionRate));  // Quantity of secretion
		trimSecretionHistory();
	}

	public Vector<SecretionRecord> getSecretionHistory()
	{	return secretionHistory;	}

	public boolean getSecretingAttractant()	{	return secretingAttractant;	}

	/** Avoid having to repeatedly perform this operation. Do it only once per time step instead, and save it. */
	public ConcentrationPerception perceiveAttractant()
	{
		if (attractantPerception == null || attractantPerception.expired())
		{
			SamplePoints points = getSamplePoints();  // Absolute-space locations around cell, to sample attractant
			attractantPerception = sampleAttractantAroundCell(points);
		}
		return attractantPerception;
	}

	/** Find concentration of attractant at the specified locations. */
	protected ConcentrationPerception sampleAttractantAroundCell(SamplePoints locations)
	{
		if (secretionMode == SecretionModeTypes.BolusCentre)
		{
			final double front = gradientFromCentre(locations.front);
			final double back = gradientFromCentre(locations.back);
			final double right = gradientFromCentre(locations.right);
			final double left = gradientFromCentre(locations.left);
			final double up = gradientFromCentre(locations.up);
			final double down = gradientFromCentre(locations.down);

			ConcentrationPerceptionImmediate cp = new ConcentrationPerceptionImmediate(
					front, back, right, left, up, down, Simulation.instance.schedule.getTime());
			return cp;
		}
		// Used primarily for testing. 
		else if (secretionMode == SecretionModeTypes.WithinBolusGradientCentre)
		{
			final double front = Simulation.space.perceiveAttractantAtPoint(locations.front);
			// Calculate vector to centre of bolus.
			final BoundedCylinder bc = ((BoundedCylinder)Simulation.instance.space);
			final Double3D spaceCentre = bc.centre;
			final Double3D vectorToCentre = spaceCentre.subtract(location);	
			// ConcentrationPerceptions are relative to cell, not in absolute space. 
			// This transformation is performed through the inverse quaternion of 'orietnation'. 
			final Double3D gradientRelativeToCell = orientation.inverse().transform(vectorToCentre);

			ConcentrationPerceptionGradientMaxSupplied cp = new ConcentrationPerceptionGradientMaxSupplied(
					front, gradientRelativeToCell, Simulation.instance.schedule.getTime());
			return cp;
		}
		else {  // Default case. 
			ConcentrationPerception cp = new ConcentrationPerceptionDeferred(Simulation.space,
					locations.front, locations.back, locations.left, locations.right, locations.up, locations.down,
					Simulation.instance.schedule.getTime());
			return cp;
		}		
	}

	private double gradientFromCentre(Double3D location)
	{
		final double beta = 1e-3;
		final double distanceToCentre = distanceFromCentre(location);
		// Exponential distribution (non-normalised)
		final double concentration = Math.exp(-beta * distanceToCentre);
		return concentration;
	}

	private double distanceFromCentre(Double3D location)
	{
		final BoundedCylinder bc = ((BoundedCylinder)Simulation.instance.space);
		final Double3D spaceCentre = bc.centre;
		final Double3D vectorToCentre = spaceCentre.subtract(location);
		final double distanceToCentre = vectorToCentre.length();
		return distanceToCentre;
	}


	/** For computational efficiency. Very old secretion histories have negligible effect in present time. */
	private void trimSecretionHistory()
	{
		if (secretionMemory == Double.POSITIVE_INFINITY)
			return;  // Never trim the history. 
		
		double time = Simulation.instance.schedule.getTime();
		Vector<SecretionRecord> delete = new Vector<SecretionRecord>();
		for (SecretionRecord h : this.secretionHistory)
			if (h.getTime() + secretionMemory < time)  // Time is expressed in minutes
				delete.add(h);
		this.secretionHistory.removeAll(delete);
	}

	public Mode followGradient()
	{
		final ConcentrationPerception chemokines = perceiveAttractant();		
		final double concentration = chemokines.maxConcentration();
				
		if (concentration >= chemokineConcentrationResponseThreshold
			&& concentration < chemokineConcentrationDesensitizationThreshold) 
		{
			return Mode.DIRECTED;
		} else {
			return Mode.UNDIRECTED;
		}
	}

	public double getDiameter()
	{ 	return diameter;	}

	public double getRadius()
	{	return diameter / 2.0;	}

	public Double3D getCurrentLocation()
	{
		if(location != null)
			return location;
		else
			return Simulation.space.getCellLocation(this);
	}
	public boolean isChemotactic()
	{
		switch(followGradient())
		{
			case DIRECTED:   return true;
			case UNDIRECTED: return false;
		}
		return false;  // Unreachable code.
	}
	
	public double getChemokineConcentrationResponseThreshold()
	{
		return chemokineConcentrationResponseThreshold;
	}

	public boolean isSecretingAttractant()
	{	return secretingAttractant;		}

	public double getAttractantSecretionRate()
	{	return attractantSecretionRate;		}

	public CellLogger.Track getLogger()
	{ 	return logger;	}

	public String getType()
	{	return "t cell";	}

	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;

		n = (Node) xPath.compile("/params/TCells/secretionMode").evaluate(params, XPathConstants.NODE);
		if (n != null)
		{
			if(n.getTextContent().equals("None")) secretionMode = SecretionModeTypes.None;
			if(n.getTextContent().equals("WithinBolus")) secretionMode = SecretionModeTypes.WithinBolus;
			if(n.getTextContent().equals("Relay")) secretionMode = SecretionModeTypes.Relay;
			if(n.getTextContent().equals("BolusCentre")) secretionMode = SecretionModeTypes.BolusCentre;
			if(n.getTextContent().equals("WithinBolusGradientCentre")) secretionMode = SecretionModeTypes.WithinBolusGradientCentre;
		}
		
		n = (Node) xPath.compile("/params/TCells/attractantSecretionRate").evaluate(params, XPathConstants.NODE);
		attractantSecretionRate = Double.parseDouble(n.getTextContent());

		n = (Node) xPath.compile("/params/TCells/relayActivationAttractionThreshold").evaluate(params, XPathConstants.NODE);
		relayActivationAttractionThreshold = Double.parseDouble(n.getTextContent());

		// Expressed in degrees, and converted to radians
		n = (Node) xPath.compile("/params/TCells/chemoSensitivityCone").evaluate(params, XPathConstants.NODE);
		perceptionCone = Double.parseDouble(n.getTextContent()) * (Math.PI / 180.);

		n = (Node) xPath.compile("/params/TCells/chemotaxisConcentrationLocation").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() > 0)  // Tag is present and set.
			chemotaxisConcentrationLocation = Double.parseDouble(n.getTextContent());

		n = (Node) xPath.compile("/params/TCells/chemotaxisConcentrationShape").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() > 0)  // Tag is present and set.
			chemotaxisConcentrationShape = Double.parseDouble(n.getTextContent());
		
		// For selecting threshold concentrations for desensitization. 
		n = (Node) xPath.compile("/params/TCells/chemotaxisConcentrationDesensitizationLocation").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() > 0)  // Tag is present and set.
			chemotaxisConcentrationDesensitizationLocation = Double.parseDouble(n.getTextContent());

		n = (Node) xPath.compile("/params/TCells/chemotaxisConcentrationDesensitizationShape").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() > 0)  // Tag is present and set.
			chemotaxisConcentrationDesensitizationShape = Double.parseDouble(n.getTextContent());
		
		
		// Parameters file expresses location as the median in absolute terms. 
		// However, LogNormalDistribution expects it as the mean of exp(distro) (which is Normally distributed).
		// Need to log(location) to convert. More convenient for user to use absolute terms.
		if (chemotaxisConcentrationLocation != 0.)
		{
			chemokineResponseDistro = new LogNormalDistribution(Math.log(chemotaxisConcentrationLocation),
																   		 chemotaxisConcentrationShape);
		} else { 
			chemokineResponseDistro = null;
		}
		if (chemotaxisConcentrationLocation != 0.)
		{
			chemokineDesensitizationDistro = new LogNormalDistribution(
					Math.log(chemotaxisConcentrationDesensitizationLocation),
					chemotaxisConcentrationDesensitizationShape);
		} else { 
			chemokineDesensitizationDistro = null;
		}
	}
}
