package soluble;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * [Largely deprecated in favour of the 2DArray version.]
 * 
 * Determines the concentration of soluble factor at some distance and time away from the secretion event. 
 * This is in standardised form - shape of the spread is Gaussian - and DOES NOT take into account the 
 * quantity secreted at the secretion event. 
 * This is to accommodate variable secretion rates, so the multiplication for secretion quantity is done
 * externally to this. 
 * 
 * As a usage, an environment tracking the location, time and quantity of every secretion event within it
 * can determine the concentration of factor at any point in time by using the distance in time and space
 * to index into this implementation. 
 * 
 * The maths to calculate concentrations based on these variables can be expensive, particularly if done
 * repeatedly. 
 * Hence, this class utilises a look up table of distances in time and space. 
 * The values of this standardised form don't change, and hence only need to be calculated once, at first
 * request. 
 *  
 * 
 * @author Mark N. Read, 2018 
 *
 */
public class HeatKernelLookupTableArrayList implements HeatKernelLookup
{

	ArrayList<ArrayList<Double>> table = new ArrayList<ArrayList<Double>>();
	final private double resolution_distance;
	final private double resolution_time;
	
	final private double diffusion_coeff;  // um2/min
	
	public HeatKernelLookupTableArrayList(double res_distance, double res_time, double diffusion_coeff) 
	{
		this.resolution_distance = res_distance;
		this.resolution_time = res_time;
		this.diffusion_coeff = diffusion_coeff;
	}
	
	/** Convert real-valued distance and time measurements into indices for the lookup table */
	private int index_distance2(double distance2)
	{	return (int) (Math.sqrt(distance2) / this.resolution_distance);	}
	
	private int index_time(double time) 
	{	return (int) (time / this.resolution_time);		}
		
	public double get(double distance2, double time)
	{
		final int i_distance2 = index_distance2(distance2);
		final int i_time = index_time(time); 
		
		if (i_distance2 >= table.size())  { 
			// No value this large has been recorded yet. 
			table.ensureCapacity(i_distance2+1);
			while (table.size() < i_distance2+1)
				table.add(null);
		}
		 
		if (table.get(i_distance2) == null)
			// This value has not been set yet. 
			table.set(i_distance2, new ArrayList<Double>());
		
		ArrayList<Double> dist = table.get(i_distance2);
		if (i_time >= dist.size())
		{
			// No value this large has been recorded yet. 
			dist.ensureCapacity(i_time+1);
			while (dist.size() < i_time+1)
				dist.add(null);
		}
		
		if (dist.get(i_time) == null)
		{
			// This value has not been set yet. Calculate then set. 
			Double value = DiffusionMaths.heatKernel3D(distance2, time, diffusion_coeff);
			dist.set(i_time, value);
		}
		
		return dist.get(i_time);
	}
}
