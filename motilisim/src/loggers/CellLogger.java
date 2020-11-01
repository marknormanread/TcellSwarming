package loggers;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.Cell;
import core.Simulation;
import movement.ChemotaxisResponsive;
import sim.engine.Schedule;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.util.Double3D;


/**
 * 
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
 * This class encapsulates the functionality of logging cell movements, and writing them to the filesystem. 
 * It contains a subclass, Track, one of which is associated with each tracked cell in the simulation. 
 * 
 * It is indented that a single simulation only ever has one instance of the CellLogger class. 
 * 
 * @author Mark N. Read
 */
public class CellLogger implements Steppable
{
	private static enum TrackType {
		VOLUME,
		SLICE,
		UNBOUNDED
	}
	private static TrackType trackMode = TrackType.VOLUME;
	public static double sliceZ = 0.;  // If in slice move, how high is the slice?
	public static double sliceDepth = 7;  // How thick is the slice?
	public static double cellLoggingInterval = 30.;  // Minutes.
	
	
	// A way of keeping hold of all the cell loggers tracking this type of cell together. 
	public static ArrayList<Track> tracks = new ArrayList<Track>();
	public static String outputDir;
		
	/**
	 * A single Track object is associated with a single Cell object. The Track is responsible for tracking 
	 * the cells's movements and compiling statistics thereof. 
	 * 
	 * @author Mark Read
	 */
	public static class Track implements Steppable
	{
		public Cell target;		
				
		public ArrayList<Double3D> positionLog = new ArrayList<Double3D>(); 
		public ArrayList<Double> timeLog = new ArrayList<Double>();  // The actual time
		public ArrayList<Long> timeIterLog = new ArrayList<Long>();  // Iterative count of timesteps so far.
		public ArrayList<Boolean> chemotacticLog = new ArrayList<Boolean>();
		
		public Track(Cell targetCell)
		{
			this.target = targetCell;	
			tracks.add(this);			
			
			Schedule sched = Simulation.instance.schedule;
			double startTime = sched.getTime();
			if (startTime < 0.0)
				startTime = 0.0;
			sched.scheduleRepeating(startTime, Simulation.loggerOrdering, this, Simulation.sampleTimeSlice_min);
		}
		
		@Override
		public void step(SimState state) 
		{	
			final double time = Simulation.instance.schedule.getTime();
			if (time >= 0.)  // stop tracking things at the epoch
			{
				positionLog.add(target.getCurrentLocation());			
				timeLog.add(Simulation.instance.schedule.getTime());
				timeIterLog.add(Simulation.instance.timeIter);
				
				if (target instanceof ChemotaxisResponsive)
				{
					chemotacticLog.add(((ChemotaxisResponsive)target).followGradient() 
										== ChemotaxisResponsive.Mode.DIRECTED);
				} else
					chemotacticLog.add(false);
			}
		}
	}	
		
	public CellLogger(String dir)
	{	this.outputDir = dir;	}
	
	/** 
	 * Allows positional data to be periodically written to FS. 
	 * For very long simulations this can facilitate checking all is well, or getting some data before cluster 
	 * compute time expires. 
	 */
	@Override
	public void step(SimState state) 
	{
		writeTrackData();
	}
	
	public void writeTrackData()
	{	
		System.out.println("Writing cell track data to the filesystem.");
		try 
		{
			BufferedWriter positionOut = setupPositionOutputFile(this.outputDir + "/_Position.csv", "Position");

			int trackID = 0;	// for writing to FS. Incremented for every encounter of a cell in imaging volume.			
			for (Track track : tracks)
			{
				/* Track ID is incremented for each new track that appears inside the imaging volume, and whenever
				 * an existing track re-enters the volume. Start it on false (as newly encountered tracks are assumed
				 * to not be in the volume).
				 */				
				boolean previouslyInsideVolume = false;
				// need a new track ID for every new track. 
				for (int t = 0; t < track.positionLog.size(); t++)
				{
					double x = track.positionLog.get(t).x;
					double y = track.positionLog.get(t).y;
					double z = track.positionLog.get(t).z;	
					final double time = track.timeLog.get(t);
					if (trackable(x, y, z, track.target))							
					{
						// check if the track ID needs to be incremented. Done for every new track, and every re-entry
						// of existing tracks into the imaging volume. 
						if (previouslyInsideVolume == false)
						{	
							trackID ++;
							previouslyInsideVolume = true;
						}	
						positionOut.write(
							x + "," + y + "," + z + "," + 
							"um," + 
							track.chemotacticLog.get(t) + "," +
							track.timeIterLog.get(t) + "," +
							trackID + "," + time + "\n");
					} else {
						previouslyInsideVolume = false;
					}
				}				
			}
			positionOut.close();
				
		} catch (IOException ex) {
			System.out.println("ERROR: exception when writing to filesystem, " + ex.toString());
		}
		System.out.println("Completed writing neutrophil track data to the filesystem.");
	}
	
	/* Determines whether the location or cell supplied should be recorded. */
	public static boolean trackable(Cell cell)
	{	final Double3D loc = cell.getCurrentLocation();
		return trackable(loc.x, loc.y, loc.z, cell);
	}
	
	/* Determines whether the location or cell supplied should be recorded. */
	public static boolean trackable(double x, double y, double z, Cell cell)
	{
		// Based purely on features of the compartment, NOT on slice or volume considerations. 
		if (! cell.recordable)  
			return false;
		
		switch (trackMode) {
			case VOLUME:				
				return Simulation.space.insideImagingVolume(x, y, z);
				
			case SLICE:
				return Simulation.space.insideImagingVolume(x, y, z) 
						&& Simulation.space.insideRecordingVolume(x, y, z, cell);
			
			case UNBOUNDED:
				return true;
				
			default:  // Unreachable code
				return false;
		}
	}
	
	private static BufferedWriter setupPositionOutputFile(String fileName, String title) throws IOException
	{
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter out = new BufferedWriter(fw);
		out.write("Position_X,Position_Y,Position_Z,Unit,Chemotactic,TimePointID,TrackID,Time\n");	
		return out;
	}
	
	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;	
		
		n = (Node) xPath.compile("/params/Recording/mode").evaluate(params, XPathConstants.NODE);
		if (n != null)
		{
			if (n.getTextContent().equals("volume"))
				trackMode = TrackType.VOLUME;
			else if (n.getTextContent().equals("slice"))
			{
				trackMode = TrackType.SLICE;
				n = (Node) xPath.compile("/params/Recording/sliceZ").evaluate(params, XPathConstants.NODE);
				sliceZ = Double.parseDouble(n.getTextContent());
				n = (Node) xPath.compile("/params/Recording/sliceDepth").evaluate(params, XPathConstants.NODE);
				sliceDepth = Double.parseDouble(n.getTextContent());
			}
			else if (n.getTextContent().equals("unbounded"))
			{
				trackMode = TrackType.UNBOUNDED;
			}
		}
	}
	
}
