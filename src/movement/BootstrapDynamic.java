package movement;

import core.MigratoryCell;
import sim.util.Double3D;
import soluble.ConcentrationPerception;
import utils.Quaternion;

/**
 *
 * Runs a combination of Bootstrap and BootstrapGradient motility depending on whether the cell
 *
 * @author Mark N. Read, 2018
 *
 */
public class BootstrapDynamic implements Orientation, Translation
{
	final Bootstrap undirected_bs;
	final BootstrapGradient directed_bs;

	private ChemotaxisResponsive.Mode currentMode;


	public BootstrapDynamic()
	{
		undirected_bs = new Bootstrap();
		directed_bs = new BootstrapGradient();
	}


	@Override
	public Double3D move(Quaternion orientation)
	{
		switch(currentMode)
		{
		case UNDIRECTED:
			return undirected_bs.move(orientation);

		case DIRECTED:
			return directed_bs.move(orientation);

		default:
			return null;
		}
	}

	@Override
	/* This is called before move */
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell)
	{
		if (cell instanceof ChemotaxisResponsive == false)
			currentMode = ChemotaxisResponsive.Mode.UNDIRECTED;
		else
			currentMode = ((ChemotaxisResponsive)cell).followGradient();  // Used in move also.

		switch(currentMode)
		{
			case UNDIRECTED:
				return undirected_bs.newOrientation(orientation, cell);

			case DIRECTED:
				return directed_bs.newOrientation(orientation, cell);

			default:
					return null;
		}
	}
}
