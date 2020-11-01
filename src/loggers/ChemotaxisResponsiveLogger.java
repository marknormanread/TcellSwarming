package loggers;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import loggers.CellLogger.Track;
import movement.ChemotaxisResponsive;
import sim.engine.SimState;
import sim.engine.Steppable;

public class ChemotaxisResponsiveLogger implements Steppable
{
	public static ChemotaxisResponsiveLogger instance = new ChemotaxisResponsiveLogger();
	
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
	
	public static void writeTrackData()
	{	
		System.out.println("Writing chemotaxis responsive threshold data to filesystem.");
		try 
		{
			BufferedWriter positionOut = setupPositionOutputFile(CellLogger.outputDir + "/_ChemotaxisResponsiveTresholds.csv", "Position");

			int trackID = 0;	// For writing to FS. Incremented for every encounter of a cell in imaging volume.			
			for (Track track : CellLogger.tracks)
			{
				/* Track ID is incremented for each new track that appears inside the imaging volume, and whenever
				 * an existing track re-enters the volume. Start it on false (as newly encountered tracks are assumed
				 * to not be in the volume).
				 */				
				boolean previouslyInsideVolume = false;
				// Need a new track ID for every new track. 
				for (int t = 0; t < track.positionLog.size(); t++)
				{
					final double x = track.positionLog.get(t).x;
					final double y = track.positionLog.get(t).y;
					final double z = track.positionLog.get(t).z;	
					if (CellLogger.trackable(x, y, z, track.target))							
					{
						// Check if the track ID needs to be incremented. 
						// Done for every new track, and every re-entry of existing tracks into the imaging volume. 
						if (previouslyInsideVolume == false)
						{	
							trackID ++;
							previouslyInsideVolume = true;
							
							positionOut.write(
								trackID 
								+ "," + 
								((ChemotaxisResponsive)track.target).getChemokineConcentrationResponseThreshold()				
								+ "\n");
						}
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
	
	private static BufferedWriter setupPositionOutputFile(String fileName, String title) throws IOException
	{
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter out = new BufferedWriter(fw);
		out.write("TrackID,chemokineConcentrationResponseThreshold\n");	
		return out;
	}
}
