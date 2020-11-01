package soluble;


/**
 * [Preamble], equivalent to the ArrayList implementation, but this is faster. 
 * The downside is that the max time/distance it can operate over must be known and set in advance.  
 * 
 * Determines the concentration of soluble factor at some distance and time away from the secretion event. 
 * This is in standardised form - shape of the spread is Gaussian - and DOES NOT take into account the 
 * quantity secreted at the secretion event. 
 * This is to accommodate variable secretion rates, so the multiplication for secretion quantity must be  
 * performed externally to this. 
 * 
 * As a usage, an environment tracking the location, time and quantity of every secretion event within it
 * can determine the concentration of factor at any point in time by using the distance in time and space
 * to index into this implementation. 
 * 
 * The maths to calculate concentrations based on these variables can be expensive, particularly if done
 * repeatedly. 
 * Hence, this class utilises a look up table of distances in time and space. 
 * The values of this standardised form don't change, and hence only need to be calculated once, at first
 * request, and can be looked up thereafter. 
 * 
 * This class discretises values for time and distance from secretion events. 
 * Upon request, given time and distance values are resolved into the bins representing this discretisation. 
 * The values calculated and stored in the bin represent the bin's centre (of distance or time), rather than
 * the value requested. 
 * This is to ensure that the transitions between binned quantities are as smooth as possible. 
 *  
 * 
 * @author Mark N. Read, 2018 
 *
 */
public class HeatKernelLookupTable2DArray implements HeatKernelLookup 
{
	double[][] table;  // Stores calculation results. [time][distance]. 
	final private double resolution_distance;  // Length of distance captured in each bin. 
	final private double resolution_time;  // Length of time captured in each bin. 
	
	final private double max_distance;  // Used to set the dimensions of the lookup table. 
	final private double max_time;
	
	final private double diffusion_coeff;  // Note the units: um2/min
	
	public HeatKernelLookupTable2DArray(double res_distance, double res_time, double diffusion_coeff,
										double max_distance, double max_time) 
	{
		this.max_distance = max_distance;
		this.max_time = max_time;
		
		int time_bins = (int) Math.ceil(max_time / res_time) + 1;
		int dist_bins = (int) Math.ceil(max_distance / res_distance) + 1;
		
		table = new double[time_bins][dist_bins];
		// Initialise table with NaN; replaced upon first lookup. 
		for (int t = 0; t < time_bins; t++)
			for (int d = 0; d < dist_bins; d++)
			{
				table[t][d] = Double.NaN;
			}
		
		this.resolution_distance = res_distance;
		this.resolution_time = res_time;
		this.diffusion_coeff = diffusion_coeff;
	}
	
	/** 
	 * Convert real-valued distance measurements into indices for the lookup table. 
	 * Supply distance squared. 
	 */
	private int index_distance2(double distance2)
	{	
		final double dist = Math.sqrt(distance2);
		assert dist < this.max_distance;
		return (int) (dist / this.resolution_distance);	
	}
	
	/** Convert real-valued distance and time measurements into indices for the lookup table */
	private int index_time(double time) 
	{	
		assert time < this.max_time;
		return (int) (time / this.resolution_time);		
	}
	
	/** Convert a time index for the lookup table into the raw time value represented by the bin centre. */
	private double time_from_index(int index)
	{	
		// Add half the bin width to place in centre of domain captured by the bin. 
		return (this.resolution_time * index) + (this.resolution_time / 2.);	
	}
	
	/** Convert lookup table distance index into raw distance squared value represented by the bin centre */
	private double distance2_from_index(int index)
	{
		// Add half the bin width to place in centre of domain captured by the bin.
		final double dist = (this.resolution_distance * index) + (this.resolution_distance / 2.); 	
		return dist * dist;
	}
	
	/** 
	 * Retrieve proportion of secreted factor perceived at some time and distance after the secretion event.
	 * This is standardised value (e.g. proportion of quantity secreted), not the absolute value. 
	 * To obtain absolute value, the returned proportion must be multiplied by the quantity that was 
	 * secreted.
	 * 
	 * Supplied distance is squared distance. 
	 */
	public double get(double distance2, double time)
	{
		// Resolve distance squared and time to a bin in the lookup table. 
		final int i_distance = index_distance2(distance2);
		final int i_time = index_time(time); 
		
		return indexed_lookup(i_distance, i_time);				
	}
	
	/**
	 * Retrieves value from the lookup table at the supplied bin (specified through indexes). 
	 * Calculates if not already stored. 
	 * @param i_distance index into lookup table. 
	 * @param i_time
	 * @return
	 */
	double indexed_lookup(int i_distance, int i_time)
	{
		try {
			final double val = table[i_time][i_distance];
	
			// Input to the heat kernel. Use values representing the centre of the bin. 
			final double bin_centred_time = time_from_index(i_time);
			final double bin_centred_distance2 = distance2_from_index(i_distance);
			// Not yet been calculated
			if (Double.isNaN(val))  
			{			
				// This value has not been set yet. Calculate then set.
				final double value = DiffusionMaths.heatKernel3D(bin_centred_distance2, bin_centred_time, 
															     diffusion_coeff);
				table[i_time][i_distance] = value;
			}
		} catch (ArrayIndexOutOfBoundsException e)
		{
			System.out.println("Requested distance and time exceed the range of values stored in the lookup"
					+ "table. This can be fixed by instantiating a larger look up table, but that needs to"
					+ "be done in code.");
		}
		
		return table[i_time][i_distance];
	}
}
