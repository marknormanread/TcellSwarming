package core;

import java.util.ArrayList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import environment.BoundedCylinder;
import loggers.CellLogger;
import loggers.ChemotaxisResponsiveLogger;
import loggers.SeedLogger;
import loggers.TimeLogger;
import sim.engine.Schedule;
import sim.util.Double3D;
import soluble.Attractant;


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
 * @author Mark N. Read, 2018
 *
 */
public class SimulationChemotacticAgent extends Simulation
{
	public static String defaultParametersPath = "parameters.xml";

	public static ArrayList<MigratoryCell> agents = new ArrayList<MigratoryCell>();

	public static int numAgents = 0;  // Number of cells in the imaging volume.
	public static int totalAgents;   // Calculated at launch. Includes cells in the buffer.

	public static enum BolusChemokineSourceType {
		NONE,
		EPOCH,  // Single secretion event at start of simulation.
		CONTINUOUS
	}
	public static BolusChemokineSourceType bolusChemokineSourceMode = BolusChemokineSourceType.NONE;
	public static double bolusChemokineSecretionRate = 0;
	
	public SimulationChemotacticAgent()
	{
		super();
		readInParameters();
	}

	public SimulationChemotacticAgent(String outDir, String paramsPath)
	{
		super(outDir, paramsPath, false);
		readInParameters();
	}
	
	public SimulationChemotacticAgent(String outDir, String paramsPath, long seed)
	{
		super(outDir, paramsPath, false, seed);
		readInParameters();
	}
	
	
	private void readInParameters()
	{
		try{
			SimulationChemotacticAgent.loadParameters(parameters);
			ChemotacticAgent.loadParameters(parameters);
			Attractant.loadParameters(parameters);
		}
		catch(XPathExpressionException e) {
			System.out.println("ERROR reading in parameters: " + e.toString());
		}		
	}
	
	public String getDefaulParametersPath()
	{	return defaultParametersPath;		}

	/**
	 * Populate the simulation's spatial environment with T cells.
	 */
	public void populateCells()
	{
		if (trackCells)
		{
			cellLogger = new CellLogger(outputPath);
			schedule.scheduleRepeating(CellLogger.cellLoggingInterval, Simulation.loggerOrdering, 
					ChemotaxisResponsiveLogger.instance, 
					CellLogger.cellLoggingInterval);
		}
		totalAgents = (int) Math.round((space.volumeSimulated()/space.volumeImaged()) * numAgents);
		System.out.println("total number of agents = " + totalAgents);

		for(int n = 0; n < totalAgents; n++)
		{
			ChemotacticAgent agent = new ChemotacticAgent(Simulation.instance.schedule);
			Simulation.space.placeCellRandomly(agent);
			agents.add(agent);
		}
		if (trackCells)
		{
			// Record neutrophil initial positions; logging only happens after neutrophils have been stepped.
			for (MigratoryCell c : agents)
				((ChemotacticAgent)c).getLogger().step(this);
		}

		// Set up a single source of chemoattractant in bolus centre, if required. 
		if (bolusChemokineSourceMode==BolusChemokineSourceType.NONE)
		{	/* Nothing needs to be done; included here for completeness. */	} 
		else {	
			
			// The following functionality relies on having a bolus, currently only implemented in BoundedCylinder.
			assert(Simulation.space instanceof BoundedCylinder); 
			assert(BoundedCylinder.bolusPresent);
			Double3D location = BoundedCylinder.centre;			
			StationaryChemokineSource chemokineSource = null;
			
			if (bolusChemokineSourceMode==BolusChemokineSourceType.EPOCH) 
			{
				chemokineSource = new StationaryChemokineSource
						(location, StationaryChemokineSource.SecretionType.EPOCH, bolusChemokineSecretionRate);
				schedule.scheduleOnce(0., Simulation.compartmentOrdering, chemokineSource);								

			}
			if (bolusChemokineSourceMode==BolusChemokineSourceType.CONTINUOUS)	
			{
				chemokineSource = new StationaryChemokineSource
						(location, StationaryChemokineSource.SecretionType.CONTINUOUS, bolusChemokineSecretionRate);
				schedule.scheduleRepeating(0., Simulation.compartmentOrdering, 
						                   chemokineSource, 
						                   Simulation.timeSlice_min);
			}
			space.cellField.setObjectLocation(chemokineSource, location);
		}
	}

	/** Tears down the simulation, can be used for writing IO. */
	@Override
	public void finish()
	{
		super.finish();
		if (trackCells)
		{
			System.out.println("Writing simulation output data to filesystem: " + outputPath);
			cellLogger.writeTrackData();			
		}
		System.out.println("Simulation completed, you may close any open windows now.");
	}


	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;
		n = (Node) xPath.compile("/params/ChemotacticAgent/numAgents").evaluate(params, XPathConstants.NODE);
		numAgents = Integer.parseInt(n.getTextContent());

		n = (Node) xPath.compile("/params/ChemotacticAgent/bolusChemokineSource").evaluate(params, XPathConstants.NODE);
		if (n.getTextContent().equals("NONE"))
			bolusChemokineSourceMode = BolusChemokineSourceType.NONE;
		if (n.getTextContent().equals("EPOCH"))
			bolusChemokineSourceMode = BolusChemokineSourceType.EPOCH;
		if (n.getTextContent().equals("CONTINUOUS")) 
			bolusChemokineSourceMode = BolusChemokineSourceType.CONTINUOUS;

		n = (Node) xPath.compile("/params/ChemotacticAgent/bolusChemokineSecretionRate").evaluate(params, XPathConstants.NODE);
		if (n != null && ! "".equals(n.getTextContent()))
			bolusChemokineSecretionRate = Double.parseDouble(n.getTextContent());		
	}


	public static void main(String[] args)
	{ 
		Simulation.cmdLineArgs = args;
		readArgs(args);
		Simulation state = new SimulationChemotacticAgent();
		execute(state);
	}
}
