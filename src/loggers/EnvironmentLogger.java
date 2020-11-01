package loggers;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import core.EpidermalCell;
import core.HairFollicle;
import core.Simulation;
import environment.BoundedCubeNeutrophil;
import sim.util.Double3D;

public class EnvironmentLogger 
{
	public static void writeEnvironmentalData(String dir)
	{
		System.out.println("Writing follicle data to the filesystem.");
		try 
		{
			BufferedWriter folliclePos = setupFolliclePositionOutputFile(dir + "/_FolliclePosition.csv");
			BufferedWriter burnPos = setupBurnPositionOutputFile(dir + "/_BurnPosition.csv");
			ArrayList<HairFollicle> follicles = ((BoundedCubeNeutrophil)Simulation.space).follicles; 
			for (int id = 0; id < follicles.size(); id++)
			{
				HairFollicle f = follicles.get(id);
				Double3D fLoc = Simulation.space.cellField.getObjectLocation(f);
				folliclePos.write(fLoc.x + "," + fLoc.y + "," + fLoc.z + "," + f.getDiameter() + ",um," + id + "\n");
			}
			
			ArrayList<EpidermalCell> burns = ((BoundedCubeNeutrophil)Simulation.instance.space).lysedEpidermals;
			for (int id = 0; id < burns.size(); id++)
			{
				EpidermalCell ec = burns.get(id);
				Double3D loc = ec.getCurrentLocation();
				burnPos.write(loc.x + "," + loc.y + "," + loc.z + ",um," + id + "\n");
			}
			folliclePos.close();
			burnPos.close();
			
			BufferedWriter simOut = setupSimulationParamsOutputFile(dir + "/_SimulationParams.csv");
			simOut.write(Simulation.timeSlice_min + "," + 
						Simulation.space.extremeHeight + "," + Simulation.space.extremeWidth + "," 
						+ Simulation.space.extremeDepth);
			
		} catch (IOException ex) {
			System.out.println("ERROR: exception when writing to filesystem, " + ex.toString());
		}
	}
	
	private static BufferedWriter setupFolliclePositionOutputFile(String fileName) throws IOException
	{
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter out = new BufferedWriter(fw);
		out.write("\n");
		out.write("Follicle Locations\n");
		out.write(" ==================== \n");
		out.write("Position X,Position Y,Position Z,Diameter,Unit,Parent\n");	
		return out;
	}
	
	private static BufferedWriter setupBurnPositionOutputFile(String fileName) throws IOException
	{
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter out = new BufferedWriter(fw);
		out.write("\n");
		out.write("Burn Locations\n");
		out.write(" ==================== \n");
		out.write("Position X,Position Y,Position Z,Unit,Parent\n");	
		return out;
	}
	
	private static BufferedWriter setupSimulationParamsOutputFile(String fileName) throws IOException
	{
		FileWriter fw = new FileWriter(fileName);
		BufferedWriter out = new BufferedWriter(fw);
		out.write("\n");
		out.write("Environmental parameters\n");
		out.write(" ==================== \n");
		out.write("timeSlice,tissueHeight,tissueWidth,tissueDepth\n");	
		return out;		
	}
}