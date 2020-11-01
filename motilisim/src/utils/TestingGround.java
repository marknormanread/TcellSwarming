package utils;

import java.util.ArrayList;
import java.util.Random;

import org.apache.commons.math3.distribution.LogNormalDistribution;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

public class TestingGround {

	public static void main(String[] args)
	{
		lookupTablePerformance();
		//logNormalDistro();
	}
	
	private static void lookupTablePerformance()
	{
		int max_times = 2280;
		// 3400 well radius + 1200 bolus radius = 4600. Div by three for resolution. 
		int max_distances = 1540;  
		ArrayList<ArrayList<Double>> al = new ArrayList<ArrayList<Double>>();
		
		double[][] ar = new double[max_times][max_distances];
		
		// Populate tables. 
		for (int t = 0; t < max_times; t++)
		{
			al.add(new ArrayList<Double>());
			ArrayList<Double> al_t = al.get(t);
			for (int d = 0; d < max_distances; d ++)
			{
				final double val = Math.random();
				ar[t][d] = val;
				al_t.add(val);			
			}
		}
		
		// Set up indexes to look up. Must be the same for each look up scheme to be a fair test. 
		int samples = (int) 1e9;
//		int[] lu_d = new int[samples];
//		int[] lu_t = new int[samples];
		Random rand = new Random();
//		for(int lu = 0; lu < samples; lu++)
//		{
//			lu_d[lu] = rand.nextInt(max_distances);
//			lu_t[lu] = rand.nextInt(max_times);
//		}
		
		// Test time for arrays. 
		System.out.println("Starting lookup for 2D array");
		// Start timer
		long startTime = System.currentTimeMillis();			
		double total = 0;
		for(int lu = 0; lu < samples; lu++)
		{
			total += ar[rand.nextInt(max_times)][rand.nextInt(max_distances)];
		}
		long endTime = System.currentTimeMillis();
		double duration = (endTime - startTime) / 1000;  // Seconds
		System.out.println("Finished; took " + duration + " seconds. total = " + total);
		
		
		// Test time for arrays. 
		System.out.println("Starting lookup for 2D ArrayLists");
		// Start timer
		startTime = System.currentTimeMillis();
		for(int lu = 0; lu < samples; lu++)
		{
			total -= al.get(rand.nextInt(max_times)).get(rand.nextInt(max_distances));
		}
		System.out.println("Finished; total = " + total);
		endTime = System.currentTimeMillis();
		duration = (endTime - startTime) / 1000;  // Seconds
		System.out.println("Finished; took " + duration + " seconds. total = " + total);
	}
	
	
	/* Checking behaviour of this distribution, cross referencing against python's scipy.stats (externally). 
	 * Need to log the location to induce median value in absolute terms (not exp(distro)). 
	 */
	private static void logNormalDistro()
	{
		
		// Parameters file expresses location as the median in absolute terms. 
		// However, LogNormalDistribution expects it as the mean of exp(distro) (which is Normally distributed).
		// Need to log(location) to convert. More convenient for user to use absolute terms. 
		LogNormalDistribution chemokineSensitivityDistro = new LogNormalDistribution(Math.log(1e-2), 4);
		
		for (int i=0; i<200; i++)
			System.out.print(chemokineSensitivityDistro.sample() + ",");	
	}
	
	private static void orientations() 
	{

		Vector3D x = new Vector3D(1, 0, 0);
		Vector3D y = new Vector3D(0, 1, 0);
		Vector3D z = new Vector3D(0, 0, 1);
		
//		Vector3D d = new Vector3D(1, 0, 0);
		
		
		/* phi = polar. Theta = azimuth */
//		SphericalCoordinates sc = new SphericalCoordinates(x);
//		System.out.println("vec = " + x + "; phi = " + sc.getPhi() + "; theta = " + sc.getTheta());
//		
//		sc = new SphericalCoordinates(y);
//		System.out.println("vec = " + y + "; phi = " + sc.getPhi() + "; theta = " + sc.getTheta());
//		
//		sc = new SphericalCoordinates(z);
//		System.out.println("vec = " + z + "; phi = " + sc.getPhi() + "; theta = " + sc.getTheta());
//		
//		sc = new SphericalCoordinates(g);
//		System.out.println("vec = " + g + "; phi = " + sc.getPhi() + "; theta = " + sc.getTheta());
		
//		SphericalCoordinates sc_ = new SphericalCoordinates(g);
		
		Vector3D g = new Vector3D(1, 1, 0);
		
		Vector3D d = g.normalize();
		
		Vector3D g__ = new Vector3D(1, 0, 0).normalize();
		
		/* The rotations defined here attempt to preserve orientation of by maintaining a notion of "up" (though any
		 * direction will do, we use the y axis). 
		 * This is needed to preserve roll autocorrelation. 
		 * Without the "up" preservation, the rotation from x axis to gradient can reorientatate y and z - there are 
		 * infinite rotations between two vectors in 3D. 
		 * The resultant reorientations of x and source gradient may not be the same as the subsequent rotation between
		 * x and the current gradient. 
		 * If this discrepancy keeps occurring, autocorrelation will be destroyed. This doesn't matter in the extreme
		 * case where a source cell follows the source gradient exactly. 
		 * There's a special case for when the gradient aligns perfectly with the vector being considered "up", as 
		 * the rotation becomes underspecified. */
		Rotation r;
		if (g.dotProduct(y) == 1)
			// vectors are parallel
			r = new Rotation(x, g);
		else
			r = new Rotation(x, y, g, y);		
		Vector3D d_ = r.applyInverseTo(d);  // d with respect to the gradient 
		System.out.println("d_ " + d_);
				
		Rotation r__;	
		double dot_prod = g__.dotProduct(y);
		if(dot_prod == 1)  // vectors are parallel
		{
			System.out.println("parallel");
			r__ = new Rotation(x, g__);
		} else { 
			r__ = new Rotation(x, y, g__, y);
		}
		Vector3D d__ = r__.applyTo(d_);
		System.out.println("g__ = " + g__);
		System.out.println("d__ = " + d__);
	}
}
