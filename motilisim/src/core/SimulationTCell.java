package core;

import java.util.ArrayList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import loggers.CellLogger;
import loggers.ChemotaxisResponsiveLogger;
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
 * @author Mark N. Read
 *
 */
public class SimulationTCell extends Simulation
{
//	public static String defaultParametersPath = "parameters-tmp.xml";
	public static String defaultParametersPath = "results/20181212-well-bootstrap_reconstruction_unbiased/parameters.xml";

	public static ArrayList<TCell> tcells = new ArrayList<TCell>();

	public static int numTCells = 50;  // Number of cells in the imaging volume.
	public static int totalTCells;   // Calculated at launch. Includes cells in the buffer.

	public SimulationTCell()
	{
		super();
		readInParameters();
	}
	
	public SimulationTCell(String outDir, String paramsPath)
	{
		super(outDir, paramsPath, false);
		readInParameters();
	}

	private void readInParameters()
	{
		try{
			SimulationTCell.loadParameters(parameters);
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
		totalTCells = (int) Math.round((space.volumeSimulated()/space.volumeImaged()) * numTCells);
		System.out.println("total number of T cells = " + totalTCells );

		for(int n = 0; n < totalTCells; n++)
		{
			TCell tcell;
			tcell = new TCell(Simulation.instance.schedule);

			Simulation.space.placeCellRandomly(tcell);
			tcells.add(tcell);
		}
		if (trackCells)
		{
			// Record initial cell positions; logging only happens after cells have been stepped.
			for (TCell c : tcells)
				c.getLogger().step(this);
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
		n = (Node) xPath.compile("/params/TCells/numTCells").evaluate(params, XPathConstants.NODE);
		numTCells = Integer.parseInt(n.getTextContent());
		TCell.loadParameters(params);
	}


	public static void main(String[] args)
	{ 
		Simulation.cmdLineArgs = args;
		readArgs(args);
		Simulation state = new SimulationTCell();
		execute(state);
	}
}
