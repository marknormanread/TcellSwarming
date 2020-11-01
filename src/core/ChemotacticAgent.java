package core;

import java.util.Vector;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.Cell.SamplePoints;
import core.MigratoryCell.OrientationParadigm;
import core.TCell.SecretionModeTypes;
import environment.BoundedCylinder;
import environment.Compartment3D.MoveResults;
import loggers.CellLogger;
import movement.ChemotaxisResponsive;
import movement.ChemotaxisResponsive.Mode;
import sim.engine.Schedule;
import sim.engine.SimState;
import sim.util.Double3D;
import soluble.AttractantSecretor;
import soluble.ConcentrationPerception;
import soluble.ConcentrationPerceptionDeferred;
import soluble.ConcentrationPerceptionGradientMaxSupplied;
import soluble.ConcentrationPerceptionImmediate;
import utils.Quaternion;

/** Motivated by Jack Hywood's swarming metric paper. 
 * This is a generic agent that can respond to chemotactic gradients in a number of manners. 
 * 
 * @author Mark N. Read
 *
 */

public class ChemotacticAgent extends MigratoryCell implements AttractantSecretor, ChemotaxisResponsive
//public class ChemotacticAgent extends TCell
{
	public static double diameter = 12.0;  // In microns.
	
	private static double attractantSecretionRate = 0.0;  // Molecules per minute.
	
	private static double chemokineGradientResponseThreshold = Double.POSITIVE_INFINITY;   
	
	// [0, 1]. Max weighting of random walk vs chemokine gradient in determining chemotaxing direction. 
	// 0 = all random walk; 1 = all chemokine gradient. 
	private static double maxChemokineAttractionWeighting;
	// Takes value >0. With large values, small chemokine gradients induce strong chemotaxing. 
	private static double heterogeneousAttractionChemokineWeighting;
	
	// Where and when the cell has been secreting
	private Vector<SecretionRecord> secretionHistory = new Vector<SecretionRecord>();
	private double secretionMemory = Double.POSITIVE_INFINITY;  // Period of recent secretion events maintained. 
	
	public static enum MotileResponseType
	{	
		UNRESPONSIVE,
		// Agent chemotactic if chemokine gradient > some threshold. 
		// Strength of chemotactic response is constant. 
		HOMOGENEOUS_ATTRACTION,   
		// Agent chemotactic if chemokine gradient > some threshold. 
		// Strength of chemotactic response related to chemokine concentration.  
		HETEROGENEOUS_ATTRACTION,
		REPULSED
	}
	protected static MotileResponseType motileResponseMode = MotileResponseType.UNRESPONSIVE;
	
	public static enum BolusResponseType
	{	
		NO_SECRETION,
		SECRETION
	}
	protected static BolusResponseType bolusResponseMode = BolusResponseType.NO_SECRETION;
	
	protected boolean secretingAttractant = false;
	public boolean motile = true;  // Arrest cells in bolus
	
	private ConcentrationPerception attractantPerception;  // Chemokine gradients and concentrations around cell
	
	private CellLogger.Track logger;
	
	public ChemotacticAgent(Schedule sched)
	{
		super(sched);
		if (Simulation.trackCells)
			logger = new CellLogger.Track(this);
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
				
				// Handle response to chemotactic gradients.
				double chemoGradient;
				switch(motileResponseMode)
				{
				case UNRESPONSIVE:
					// No adjustment to orientation based on chemokines.
					break;
				case HOMOGENEOUS_ATTRACTION:														
					// Process influence of chemokine, if concentration above threshold.
					chemoGradient = perceiveAttractant().maxConcentrationDifferential();
					if (chemoGradient >= chemokineGradientResponseThreshold)
					{
						chemotacticReorientation(maxChemokineAttractionWeighting, false);
					}
					// Else, orientation is unchanged from above. 					
					break;
				case HETEROGENEOUS_ATTRACTION:
					// Process influence of chemokine, if concentration above threshold.
					chemoGradient = perceiveAttractant().maxConcentrationDifferential();
					if (chemoGradient >= chemokineGradientResponseThreshold &&
						/* Avoid log(0), in case chemoThreshold <= 0  */
						chemoGradient >= 0. )
					{
						double chemoInfluence = maxChemokineAttractionWeighting * 
								(1 - Math.exp(-heterogeneousAttractionChemokineWeighting * chemoGradient));
						chemotacticReorientation(chemoInfluence, false);
					}
					// Else, orientation is unchanged from above.
					break;
				case REPULSED:
					// Process influence of chemokine, if concentration above threshold.
					chemoGradient = perceiveAttractant().maxConcentrationDifferential();
					if (chemoGradient >= chemokineGradientResponseThreshold &&
						/* Avoid log(0), in case chemoThreshold <= 0  */
						chemoGradient > 0. )
					{
						double chemoInfluence = maxChemokineAttractionWeighting * 
								(1 - Math.exp(-heterogeneousAttractionChemokineWeighting * chemoGradient));
						chemotacticReorientation(chemoInfluence, true);
						
						
						// TODO DEBUG
//						System.out.println("CT: " + Simulation.instance.schedule.getTime() + 
//								" dist from bolus: " + distanceFromCentre(this.location) + 
//								" chemoInfluence: " + chemoInfluence);
					}
					// Else, orientation is unchanged from above. 	
					break;
				}
			}
			Double3D move = null;
			move = translationActuator.move(orientation);
			
			// Performs both intercellular and boundary collision detection.
			MoveResults mr = Simulation.space.moveCellCollisionDetection(this, move);
			location = mr.newLocation;
			bounce = mr.bounce;	 // Bounce off other cells that may have been contacted.			
		}

		// If in bolus start secreting.
		if (Simulation.space.insideBolus(location) && bolusResponseMode == BolusResponseType.SECRETION)			
			secretingAttractant = true;

		if (secretingAttractant)
			secrete((Simulation)state);
	}
	
	/**
	 * Re-orients the cell based on the chemokine gradient.  
	 * @param chemoInfluence: Extent to which chemokine gradient influences direction; 0 <= x <= 1
	 * 		 newDir = chemoInfluence.chemoDir + (1-chemoInfluence).rwDir = unit vector.  			 
	 *		 If 0.5, chemo and random walk directions have equal contribution. 
	 *		 If 0, chemokine has no influence. 
	 *		 If 1, random walk has no influence. 
	 */
	private void chemotacticReorientation(double chemoInfluence, boolean repulse)
	{
		// Process influence of chemokine, if there is a non-zero gradient. 
		// Direction of gradient at cell's current location, relative to the cell.
		Double3D chemoDir = perceiveAttractant().gradientDirection();
		if (repulse)
			chemoDir = new Double3D(-chemoDir.x, -chemoDir.y, -chemoDir.z);
		
		if (chemoDir.lengthSq() > 0.)  // Avoid div by zero if zero length (no gradient).
		{
			// Get current orientation as a unit vector in absolute space.					
			Double3D rwDir = orientation.transform(x_axis);  // x_axis is already a unit vector.
			
			chemoDir = chemoDir.normalize();  // Turn into a unit vector;
			// Transform to a gradient direction in absolute space
			chemoDir = orientation.transform(chemoDir);  
		
			/* newDir = a.chemoDir + (1-a).rwDir = unit vector. 
			 * If 0.5, chemo and random walk directions have equal contribution. 
			 * If 0, chemokine has no influence. 
			 * If 1, random walk has no influence. */ 
			chemoDir = chemoDir.multiply(chemoInfluence);
			rwDir = rwDir.multiply(1. - chemoInfluence);
			
			Double3D newDir = rwDir.add(chemoDir);
			orientation = Quaternion.faceVector(newDir);
		} 
		// Else, orientation is unchanged.
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
		ConcentrationPerception cp = new ConcentrationPerceptionDeferred(Simulation.space,
				locations.front, locations.back, locations.left, locations.right, locations.up, locations.down,
				Simulation.instance.schedule.getTime());

		return cp;
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
				
		if (concentration >= chemokineGradientResponseThreshold) {
			return Mode.DIRECTED;
		} else {
			return Mode.UNDIRECTED;
		}
	}
	
	private double distanceFromCentre(Double3D location)
	{
		final Double3D spaceCentre = BoundedCylinder.centre;
		final Double3D vectorToCentre = spaceCentre.subtract(location);
		final double distanceToCentre = vectorToCentre.length();
		return distanceToCentre;
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
		boolean result = false;
		switch(followGradient())
		{
			case DIRECTED:   
				result = true; 
				break;
			case UNDIRECTED: 
				result = false;
				break;
		}
		return result;
	}	
	
	@Override
	public double getChemokineConcentrationResponseThreshold()
	{	return chemokineGradientResponseThreshold;		}

	@Override
	public boolean isSecretingAttractant()
	{	return secretingAttractant;		}

	@Override
	public double getAttractantSecretionRate()
	{	return attractantSecretionRate;		}

	public CellLogger.Track getLogger()
	{ 	return logger;	}
	
	@Override
	public String getType()
	{	return "chemotactic agent";	}
	
	/** Load in values from parameters file pertinent to this class. */
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;

		n = (Node) xPath.compile("/params/ChemotacticAgent/motileResponseMode").evaluate(params, XPathConstants.NODE);
		if (n != null)
		{
			if(n.getTextContent().equals("UNRESPONSIVE"))
				motileResponseMode = MotileResponseType.UNRESPONSIVE;
			if(n.getTextContent().equals("HOMOGENEOUS_ATTRACTION")) 
				motileResponseMode = MotileResponseType.HOMOGENEOUS_ATTRACTION;
			if(n.getTextContent().equals("HETEROGENEOUS_ATTRACTION")) 
				motileResponseMode = MotileResponseType.HETEROGENEOUS_ATTRACTION;
			if(n.getTextContent().equals("REPULSED")) 
				motileResponseMode = MotileResponseType.REPULSED;			
		}
		
		n = (Node) xPath.compile("/params/ChemotacticAgent/bolusResponseMode").evaluate(params, XPathConstants.NODE);
		if (n != null)
		{
			if(n.getTextContent().equals("NO_SECRETION"))
				bolusResponseMode = BolusResponseType.NO_SECRETION;
			if(n.getTextContent().equals("SECRETION")) 
				bolusResponseMode = BolusResponseType.SECRETION;			
		}		
		
		if (bolusResponseMode == BolusResponseType.SECRETION)
		{
			n = (Node) xPath.compile("/params/ChemotacticAgent/attractantSecretionRate").evaluate(params, XPathConstants.NODE);
			attractantSecretionRate = Double.parseDouble(n.getTextContent());
		}
		
		// Expressed in degrees, and converted to radians
		n = (Node) xPath.compile("/params/ChemotacticAgent/chemokineConcentrationResponseThreshold").evaluate(params, XPathConstants.NODE);
		chemokineGradientResponseThreshold = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/ChemotacticAgent/maxChemokineAttractionWeighting").evaluate(params, XPathConstants.NODE);
		maxChemokineAttractionWeighting = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/ChemotacticAgent/heterogeneousAttractionChemokineWeighting").evaluate(params, XPathConstants.NODE);		
		heterogeneousAttractionChemokineWeighting = Double.parseDouble(n.getTextContent());
	}
}
