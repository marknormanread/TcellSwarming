package loggers;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import core.Simulation;

/** 
 * Record the random number generator seed, for reproduciblity. 
 * 
 * @author Mark N. Read, 2018
 *
 */
public class SeedLogger 
{
	public static void writeSeedData(String dir)
	{
		System.out.println("Writing random number seed to the filesystem.");
		try
		{
			String fileName = dir + "/_Seed.txt";
			FileWriter fw = new FileWriter(fileName);
			BufferedWriter timeOut = new BufferedWriter(fw);
			// convert from minutes to seconds
			timeOut.write(Long.toString(Simulation.seed));			
			timeOut.close();
				
		} catch (IOException ex) {
			System.out.println("ERROR: exception when writing to filesystem," + ex.toString());
		}		
	}
}
