package core;

import java.util.ArrayList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import loggers.CellLogger;
import loggers.EnvironmentLogger;
import loggers.SwarmLogger;
import loggers.TimeLogger;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import environment.BoundedCubeNeutrophil;
import environment.Compartment3D;
import soluble.Attractant;
import soluble.LTB4;

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
    
 * Entry point to a simulation of neutrophil motility. 
 * 
 * 
 * @author Mark N. Read
 *
 */
public class SimulationNeutrophil extends Simulation 
{
	public static String defaultParametersPath = "parameters.xml";
	
	public static ArrayList<Neutrophil> neutrophils = new ArrayList<Neutrophil>();	
	
	public static SwarmLogger swarmLogger = null;
	
	public static int numNeutrophils = 50;  // Number of neutrophils in the imaging volume.
	public static int totalNeutrophils;  // Calculated at launch. Includes neutrophils in the buffer.	
	public static double follicleRatio = 0.1;
	
	public static boolean burn = true;
	public static double burnTime = 0.0;
	
	public SimulationNeutrophil()
	{	
		super();
		readInParameters();
    }
	
	public SimulationNeutrophil(String outDir, String paramsPath)
	{	
		super(outDir, paramsPath, false);
		readInParameters();
    }
		
	private void readInParameters() 
	{
		try{
			SimulationNeutrophil.loadParameters(parameters);
			Neutrophil.loadParameters(parameters);
			LTB4.loadParameters(parameters);
			Attractant.loadParameters(parameters);
			SwarmLogger.loadParameters(parameters);
		}
		catch(XPathExpressionException e) { 
			System.out.println("ERROR reading in parameters: " + e.toString());
		}
	}
	/**
	 * Subclasses can override this if a different compartment is needed. 
	 */
	public Compartment3D initializeCompartment()
	{	
		return new BoundedCubeNeutrophil(true, false);		
	}
	
	public String getDefaulParametersPath()
	{	return defaultParametersPath;		}
	
	/**
	 * Populate the simulation's spatial environment with neutrophils.
	 */
	public void populateCells()
	{	
		// create a new laser blast, which schedules itself for the given time. 
		if (burn)
			schedule.scheduleOnce(burnTime, blastOrdering, new LaserBlast(burnTime));
		
		if (trackCells)
		{
			cellLogger = new CellLogger(outputPath);
			swarmLogger = new SwarmLogger();
			schedule.scheduleRepeating(0.0, Simulation.loggerOrdering, swarmLogger, timeSlice_min);
		}
			
		totalNeutrophils = (int) Math.round( (space.volumeSimulated()/space.volumeImaged()) * numNeutrophils);
		System.out.println("total number of neutrophils simulated = " + totalNeutrophils );
		
		for(int n = 0; n < totalNeutrophils; n++)
		{
			Neutrophil neutro = new Neutrophil(Simulation.instance.schedule);
			Simulation.space.placeCellRandomly(neutro);
			
			neutrophils.add(neutro);
		}
		if (trackCells)
		{
			// Record neutrophil initial positions; logging only happens after neutrophils have been stepped.
			for (Neutrophil n : neutrophils)
				n.getLogger().step(this);
		}
	}
	
	/** Tears down the simulation, can be used for writing IO. */ 
	public void finish()
	{
		super.finish();
		if (trackCells)
		{
			System.out.println("Writing simulation output data to filesystem: " + outputPath);
			cellLogger.writeTrackData();
			swarmLogger.writeSwarmData(outputPath);
			EnvironmentLogger.writeEnvironmentalData(outputPath);
		}
		System.out.println("Simulation completed, you may close any open windows now.");
	}

	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;
		n = (Node) xPath.compile("/params/Neutrophil/numNeutrophils").evaluate(params, XPathConstants.NODE);
		numNeutrophils = Integer.parseInt(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Neutrophil/burn")
				.evaluate(params, XPathConstants.NODE);		
		burn = Boolean.parseBoolean(n.getTextContent());
		n = (Node) xPath.compile("/params/Neutrophil/burnTime")
				.evaluate(params, XPathConstants.NODE);		
		burnTime = Double.parseDouble(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Neutrophil/follicleRatio")
				.evaluate(params, XPathConstants.NODE);		
		follicleRatio = Double.parseDouble(n.getTextContent());
	}
	
	
	public static void main(String[] args)
	{
		readArgs(args);
		Simulation state = new SimulationNeutrophil();
		execute(state);
	}
}
