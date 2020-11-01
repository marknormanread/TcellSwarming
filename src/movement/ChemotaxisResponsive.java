package movement;

/** 
 * Implemented by cells that respond to chemotactic gradients. 
 * 
 * @author Mark N. Read, 2018
 *
 */
public interface ChemotaxisResponsive 
{
	public static enum Mode
	{
		DIRECTED,
		UNDIRECTED,		
	};
	
	public Mode followGradient();	
	
	public double getChemokineConcentrationResponseThreshold();
}
