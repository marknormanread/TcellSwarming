package soluble;


/**
 * Determines the concentration of chemokine at a given distance from a single secretor, a given duration 
 * after it started secreting.  
 * Inherent herein is the assumption that the secretor is stationary.
 * Another assumption is that the secretor is still secreting at the current point in time. 
 * Together these assumptions permit use of a lookup table that descretises the domain of this function, and
 * as such a given calculation need be performed once only.
 * 
 * A third assumption (unrelated to computational efficiency): chemokine does not decay. 
 * To be implemented, each heat-kernel calculation (a given time and distance from seceretor) would need to
 * be decayed individually, as the magnitude of the decay is time-dependent. 
 * 
 * Implementation-wise, this makes use of another lookup table for heat kernel calcuations for individual
 * time points. 
 * All this class does is record the cumulative effect of repeated secretions over a duration of time. 
 * 
 * @author Mark N. Read, 2018
 *
 */
public class StationaryDiffusionTimeIntegration implements HeatKernelLookup
{
	
	double[][] cumulativeTable;  // Look up table. [duration][distance]
	final private double resolution_distance;
	final private double resolution_duration;
	
	final private double max_distance;
	final private double max_duration;
		
	/* Stores non-cumulative standardised (secretion rate must be accounted for externally) values for 
	 * concentration at given distance and time point after a SINGLE secretion event. 
	 * This present class performs an integration over time to this data. */
	final HeatKernelLookupTable2DArray instantaneousLookup;
	
	public StationaryDiffusionTimeIntegration(double res_distance, double res_duration, 
			double diffusion_coeff, double max_distance, double max_duration)
	{
		this.max_distance = max_distance;
		this.max_duration = max_duration;
		
		int duration_bins = (int) Math.ceil(max_duration / res_duration) + 1;
		int dist_bins     = (int) Math.ceil(max_distance / res_distance) + 1;
		
		cumulativeTable = new double[duration_bins][dist_bins];
		// Initialise table with NaN; replaced upon first lookup. 
		for (int t = 0; t < duration_bins; t++)
			for (int d = 0; d < dist_bins; d++)
			{
				cumulativeTable[t][d] = Double.NaN;
			}
		
		this.resolution_distance = res_distance;
		this.resolution_duration = res_duration;
		
		// Ensure similar dimensions for cumulative and instantaneous heat kernel calculation tables. 
		instantaneousLookup = new HeatKernelLookupTable2DArray(res_distance, res_duration, diffusion_coeff,
															   max_distance, max_duration);
	}
	
	/** Convert real-valued distance squared measurements into indices for the lookup table */
	private int index_distance2(double distance2)
	{	
		// Indexing should be done on standard distance, else bin widths vary with distance, with fine 
		// resolution at large distances and poor resolution at small distances. 
		final double dist = Math.sqrt(distance2);
		assert dist < this.max_distance;  // Safety. 
		return (int) (dist / this.resolution_distance);	
	}
	
	/** Convert real-valued secretion duration measurements into indices for the lookup table */
	private int index_duration(double duration) 
	{	
		assert duration < this.max_duration;  // Safety. 
		return (int) (duration / this.resolution_duration);		
	}
	
	/**
	 * time represents how long the secretor has been secreting for. 
	 * distance2 is distance from the secretor squared. 
	 */
	public double get(double distance2, double duration)
	{
		final int i_distance = index_distance2(distance2);
		final int i_duration = index_duration(duration);
		
		return recursiveCalculator(i_distance, i_duration);		
	}
	
	/**
	 * The core of this class. 
	 * Given this class's recording of concentrations arising from a single secretor at a given duration and 
	 * distance away as cumulative, the heat kernel calculations can be iteratively updated. 
	 * Assume a given distance away from the secretor.'
	 * For duration t+1, the concentration perceived is that of duration t, plus a single additional heat
	 * kernel calculation. 
	 * Hence, for a given distance, this method recurses through iteratively shorter durations until it finds
	 * either 
	 * 1. A calculation that has been performed. In this case, it then completes and stores the calculations  
	 *    up to the requested duration.  
	 * 2. hits duration 0, at which point it likewise completes and stores calculations. 
	 * 
	 * @param i_distance Indexes into a bin in the cumulative look up table. 
	 * @param i_duration
	 * @return
	 */
	private double recursiveCalculator(int i_distance, int i_duration)
	{
		// Adds a single Heat Kernel calculation to existing calculations. 
		double val = cumulativeTable[i_duration][i_distance];
		
		// Not yet been calculated
		if (Double.isNaN(val))  
		{
			val = instantaneousLookup.indexed_lookup(i_distance, i_duration);
			if (i_duration > 0)  // Check for hitting the end of the look up table. 
				val += recursiveCalculator(i_distance, i_duration-1);
		
			cumulativeTable[i_duration][i_distance] = val;	
		}
		
		return val;
	}
	
}
