package core;


import java.util.Vector;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import environment.BoundedCubeNeutrophil;
import environment.Compartment3D.MoveResults;
import loggers.CellLogger;
import sim.engine.Schedule;
import sim.engine.SimState;
import sim.util.Double3D;
import soluble.ConcentrationPerception;
import soluble.ConcentrationPerceptionDeferred;
import soluble.ConcentrationPerceptionImmediate;
import soluble.LTB4Secretor;
import utils.Utils;

@SuppressWarnings("serial")
public class Neutrophil extends MigratoryCell implements LTB4Secretor 
{	
	public static double diameter = 13.0;  // in micrometers. 

	// The quantity of attractant required to bring a neutrophil from random walk state to attractant_stimulated. 
	private static double attractantStimulated_AttractantThreshold = Double.NaN;  // 0.1
	
	// The quantity of attractant required to bring a neutrophil from the ATTRACTANT_STIMULATED to 
	// LTB4_STIMULATED state. 
	private static double ltb4Stimualted_AttractantThreshold = Double.NaN; // 10
	// Threshold quantity of LTB4 required to bring a neutrophil from the ATTRACTANT_STIMULATED to 
	// LTB4_STIMULATED state.  
	private static double ltb4Stimualted_LTB4Threshold = Double.NaN;

	private static double ltb4SecretionRate = Double.NaN;  // Per minute  10000.0
	
	// How far off the x-axis a cell can perceive a gradient. In PI-radians (1 PI-radian = 180 degrees).  
	private static double perceptionCone = 0.33333 * Math.PI;
	
	private ConcentrationPerception attractantPerception;
	private ConcentrationPerception ltb4Perception;
	
	private Vector<SecretionRecord> secretionHistory = new Vector<SecretionRecord>();
		
	private CellLogger.Track logger;
	// States the cell can exist in. 
	private NeutrophilState cellState;
	private enum NeutrophilState {
		RANDOM_WALK, ATTRACTANT_STIMUALTED, LTB4_STIMUALTED
	}
	
	private boolean inSwarm = false; 
	 
	/* Don't use this, only here to make RenderableNeutrophil work. */
	public Neutrophil()
	{}
	
	public Neutrophil(Schedule sched)
	{
		super(sched);				
		cellState = NeutrophilState.RANDOM_WALK;
		if (Simulation.trackCells)
			logger = new CellLogger.Track(this);
	}
	
	@Override
	public void step(SimState state) 
	{	// hello i am going to make changes 
		location = Simulation.space.getCellLocation(this); 
		if (bounce.lengthSq() == 0.0)
		{
			attractantPerception = perceiveAttractant();			
			ltb4Perception = perceiveLTB4();
			if (cellState == NeutrophilState.RANDOM_WALK)
			{	
				// Sufficient LTB4 to become LTB4-stimulated?
				if (ltb4Perception.maxConcentration() > ltb4Stimualted_LTB4Threshold)
					cellState = NeutrophilState.LTB4_STIMUALTED;
				// If not, sufficient attractant to become attractant-stimulated?
				else if (attractantPerception.maxConcentration() > attractantStimulated_AttractantThreshold)
					cellState = NeutrophilState.ATTRACTANT_STIMUALTED;
				else 
					// No change in neutrophil state; perform random walk.
					orientation = orientationActuator.newOrientation(orientation, this);						
			}
			
			if (cellState == NeutrophilState.ATTRACTANT_STIMUALTED)
			{
				if (attractantPerception.maxConcentration() > ltb4Stimualted_AttractantThreshold)
					cellState = NeutrophilState.LTB4_STIMUALTED;				
				else  // No change in neutrophil state. 
				{
					rotateToVectorLimited(attractantPerception.gradientDirection(), perceptionCone);
				}
			}
			
			if (cellState == NeutrophilState.LTB4_STIMUALTED)
			{
				// Convert to unit vectors to ensure that each gradient has equal influence on the cell. 
				Double3D attractantGradient = Utils.unitVector(attractantPerception.gradientDirection());
				Double3D ltb4Gradient = Utils.unitVector(ltb4Perception.gradientDirection());
				Double3D direction = attractantGradient.add(ltb4Gradient);
				rotateToVectorLimited(direction, perceptionCone);
				
				double secretionRate = ltb4SecretionRate * Simulation.timeSlice_min;
				// ----- TEST -------
				// DOES VARIABLE SECRETION RATE IMPROVE DYNAMICS?
	////			if (attractantPerception.getMax() > ltb4Stimualted_AttractantThreshold)
	////				secretionRate = ltb4SecretionRate * Simulation.timeslice;
	////			else
	////			{
	////				double mean = ltb4Perception.getMean();
	////				double prop = Math.log(mean + 1.0) / Math.log(ltb4SecretionRate+1);
	////				double compulsory = 0.0;//ltb4SecretionRate * 0.01;
	////				double variable = ltb4SecretionRate - compulsory;
	////				secretionRate = (prop * variable) + compulsory;
	////				secretionRate *= Simulation.timeslice;
	////			}
				// ----- END TEST -----
				
				// Get location of back of cell (negative radius along x-axis) and transform to absolute space.  
				final Double3D back = orientation.transform(-getRadius(), 0.0, 0.0).add(location);
				secretionHistory.add(new SecretionRecord(
						back,							// location where secretion takes place. Back of cell. 
						state.schedule.getTime(),		// time of secretion
						secretionRate));				// quantity of secretion
				
				trimSecretionHistory();
			}				
		}
		bounce();
		Double3D move = null;
		move = translationActuator.move(orientation);		
		
		// Used for debugging purposes.
		//move = calculateMoveConstant();
		MoveResults mr = Simulation.space.moveCellCollisionDetection(this, move);
		location = mr.newLocation;
		bounce = mr.bounce;	 // Bounce off other cells that may have been contacted.
		collidedCells = mr.colliders;
	}
	
	/* Used for debugging, has cell move at a constant speed. */
	private Double3D calculateMoveConstant()
	{
		Double3D facing = orientation.transform(x_axis);
		double dist = 5.0;
		return new Double3D(facing.x*dist, facing.y*dist, facing.z*dist);
	}
			
	/** Avoid having to repeatedly perform this operation. Do it only once per time step instead, and save it. */
	public ConcentrationPerception perceiveAttractant()
	{
		if (attractantPerception == null || attractantPerception.expired())
		{
			SamplePoints points = getSamplePoints();  // Absolute-space locations around cell, to sample attractant
			attractantPerception = new ConcentrationPerceptionDeferred(Simulation.space, 
					points.front, points.back, points.left, points.right, points.up, points.down,
					Simulation.instance.schedule.getTime());
		}
		return attractantPerception;
	}
	
	public ConcentrationPerception perceiveLTB4()
	{
		if (ltb4Perception == null || ltb4Perception.expired())
		{
			SamplePoints points = getSamplePoints();  // Absolute-space locations around cell, sample soluble factor
			BoundedCubeNeutrophil space = (BoundedCubeNeutrophil)Simulation.instance.space;
			final double front = space.perceiveLTB4AtPoint(points.front);
			final double back  = space.perceiveLTB4AtPoint(points.back);
			final double right = space.perceiveLTB4AtPoint(points.right);
			final double left  = space.perceiveLTB4AtPoint(points.left);
			final double up    = space.perceiveLTB4AtPoint(points.up);
			final double down  = space.perceiveLTB4AtPoint(points.down);
			
			ltb4Perception = new ConcentrationPerceptionImmediate(
					front, back, right, left, up, down, Simulation.instance.schedule.getTime());
		}
		return ltb4Perception;
	}
		
	private void trimSecretionHistory()
	{
		double time = Simulation.instance.schedule.getTime();
		Vector<SecretionRecord> delete = new Vector<SecretionRecord>();
		for (SecretionRecord h : this.secretionHistory)
			if (h.getTime() + 50.0 < time)
				delete.add(h);
		this.secretionHistory.removeAll(delete);
	}
	
	public double getDiameter()
	{ 	return diameter;	}
	
	public double getRadius()
	{	return diameter / 2.0;	}
	
	public CellLogger.Track getLogger()
	{ 	return logger;	}
	
	public boolean isSecretingLTB4()
	{	return cellState == NeutrophilState.LTB4_STIMUALTED;	}
		
	public double getLTB4SecretionRate()
	{	return ltb4SecretionRate;	}
	
	public Vector<SecretionRecord> getSecretionHistory()
	{	return secretionHistory;	}
	
	public static double getMaxLTB4SecretionRate()
	{	return ltb4SecretionRate;	}
	
	public boolean isStateRandomWalk()
	{	return cellState == NeutrophilState.RANDOM_WALK;	}
	
	public boolean isStateAttractantStimulated()
	{ 	return cellState == NeutrophilState.ATTRACTANT_STIMUALTED;	}
	
	public boolean isStateLTB4Stimulated()
	{	return cellState == NeutrophilState.LTB4_STIMUALTED;	}
	
	public boolean isInSwarm()
	{	return inSwarm; 	}
	
	public void setInSwarm(boolean inSwarm)
	{	this.inSwarm = inSwarm;	}
	
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;
		n = (Node) xPath.compile("/params/Neutrophil/diameter")
				.evaluate(params, XPathConstants.NODE);
		diameter = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Neutrophil/attractantStimulated_AttractantThreshold")
				.evaluate(params, XPathConstants.NODE);
		attractantStimulated_AttractantThreshold = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Neutrophil/ltb4Stimualted_AttractantThreshold")
				.evaluate(params, XPathConstants.NODE);
		ltb4Stimualted_AttractantThreshold = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Neutrophil/ltb4Stimualted_LTB4Threshold")
				.evaluate(params, XPathConstants.NODE);
		ltb4Stimualted_LTB4Threshold = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Neutrophil/ltb4SecretionRate")
				.evaluate(params, XPathConstants.NODE);
		ltb4SecretionRate = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Neutrophil/perceptionCone")
				.evaluate(params, XPathConstants.NODE);
		perceptionCone = Double.parseDouble(n.getTextContent());
	}

	public String getType()
	{	return "neutrophil";	}
}
