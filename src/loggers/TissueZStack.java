package loggers;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

import sim.engine.SimState;
import sim.engine.Steppable;
import sim.util.Double3D;
import core.Simulation;
import environment.Compartment3D;


/**
 * Class takes a picture of tissue quality, through a 2D slice of the environment. Intended to show chemokine 
 * concentrations.  
 * @author Mark N. Read, 2017
 *
 */
public abstract class TissueZStack implements Steppable 
{
	private static Color cellColour = new Color(0.0f, 0.8f, 0.0f, 1.0f);
	private Color bolusColor = new Color(0.2f, 0.2f, 0.2f, 1.0f);
	protected String dir;		// directory where files are to be written. 
	/*
	 * Abstract methods that concrete subclasses use to implement imaging of specific simulation factors. 
	 */
	abstract double getConcentration(Double3D location);
	abstract File getImageFile();
	
	/** Called by the schedule to take a snap shot. */
	@Override
	public void step(SimState state) 
	{
		Compartment3D comp = (Compartment3D)Simulation.space;
		final int resolution = 100;
		final double z = Simulation.space.extremeDepth / 2.;  // tissue depth at which snapshot is to be taken.		
		// given tissue dimension and resolution, calculate how many um each pixel represents. 
		double step = Simulation.space.extremeWidth / (float)resolution;
		BufferedImage img = new BufferedImage(resolution, resolution, BufferedImage.TYPE_INT_RGB);
 
		double[][] concentrations = new double[resolution][resolution];
		double maxConcentration = 0.0;
		// retrieve all concentrations, and find the maximum value such that full dynamic range on images is used.		
		for (int x=0; x<resolution; x++)
			for (int y=0; y<resolution; y++)
			{
				double lx = x * step;
				double ly = y * step;
				Double3D pixelLoc = new Double3D(lx, ly, z);
				concentrations[x][y] = getConcentration(pixelLoc);
				if (concentrations[x][y] > maxConcentration)
					maxConcentration = concentrations[x][y];
			}
		  
		// compile the image. 
		for (int x=0; x<resolution; x++)
			for (int y=0; y<resolution; y++)
			{
				double lx = x * step;
				double ly = y * step;
				final Double3D pixelLoc = new Double3D(lx, ly, z);				
				 
				// if pixel is within bolus
				if(comp.insideCell(pixelLoc)) {  
					img.setRGB(x, y, cellColour.getRGB());
					// pixel is within cell
				} else if (comp.insideBolus(pixelLoc)) {  
					img.setRGB(x, y, bolusColor.getRGB());					
				// could be open space, in which case show gradient.
				} else {
					if (concentrations[x][y] > 0.)
					{
						int col = convertConcentrationToColor(concentrations[x][y], maxConcentration);
						img.setRGB(x, y, col);
					}					
				}
			}
		File f = getImageFile();
		try {
			ImageIO.write(img, "PNG", f);
		} catch (IOException e) 
		{ 	e.printStackTrace(); 	}
	}
	
	/* Converts the supplied `concentration` into a greyscale color, with pure white set to occur at a value of `max`.
	 * Implements a log-scale, hence adding dynamic range to smaller values at the expensve of larger values.  
	 */
	private int convertConcentrationToColor(double concentration, double max)
	{
		double prop = concentration / max;
		prop = Math.sin(prop * Math.PI/2.);  // log-like, but maintains a range of [0, 1]
		prop = Math.sin(prop * Math.PI/2.);  // do it twice to really emphasize small quantities 
		// used for debugging
		//System.out.println("max = " + max + "; concentration = " + concentration + "; prop = " + prop);
		
		int iprop = (int) (prop * 255.0);  
		int r = iprop;	// RGB colors.  
		int g = iprop;
		int b = iprop;
		// perform shifting operations to compile the color. 
		int col = (r << 16) | (g << 8) | b;
		return col;
	}
}