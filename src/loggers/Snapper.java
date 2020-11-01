package loggers;

import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import sim.display3d.Display3D;
import sim.engine.SimState;
import sim.engine.Steppable;
import sim.util.media.PNGEncoder;
import core.Simulation;


/**
 * This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * This will create a series of png snapshots of the tissue volume in the experiment directory, under "movie".
 *
 * images are named "img_xxx.png"
 *
 * to create a movie of this, use ffmpeg (a bash utility). Run this command from the "movie" directory.
 *
 * $> ffmpeg -framerate 5 -i img_%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
 *
 *
 * @author Mark N. Read
 *
 */

public class Snapper implements Steppable
{
	private Display3D display;
	private String movieDir;
	private int sequenceNum = 0;		// used to label the snapshots taken of imaging volume sequentially.
	
	/* There's a problem with simulation timing when reenacting. 
	 * Reenacted snapped images are recorded as being one timestep later than they should be (and are when 
	 * visualising directly in a normal non-reenacted simulation). 
	 * I'm not sure what the problem is, but I think it relates to the ordering in which cells are stepped, the 
	 * simulation time is stepped, the GUI is drawn, and cell locations are drawn. 
	 * It's like the GUI is updated, and then the simulation time is stepped in antipication of the next action 
	 * (but not drawn on the GUI).
	 * I haven't been able to correct this behaviour. 
	 * This is a hack, but it works.
	 * We set the snapper to record times as one timestep duration less than the simulation's schedule reports for
	 * reenactments. 
	 */
	public static boolean reenacting = false;
	
	public Snapper(Display3D disp)
	{
		this.display = disp;
		movieDir = Simulation.outputPath + "/stills";
		File dir = new File(movieDir);
		if (!dir.exists())
		{
			System.out.println("creating directory : " + movieDir);
			boolean result = dir.mkdirs();
			if (!result)	System.out.println("ERROR: could not create directory " + movieDir);
		}
	}


	@Override
	public void step(SimState state)
	{
		/* A complete hack, and an annoying one. There seems to be some synchronization issues between the Display3D's
		 * rendering of the tissue space, and the simulation itself. This sleep is purely to allow some catch up time.
		 * Without it, some of the cells in the simulation don't move between subsequent png files generated, it looks
		 * like some cells are being moved whilst the image is still being captured. I can't find a clean way to
		 * fix this, and have spent enough time trying. Hence, sleep for a while.
		 *
		 * Annoyingly, MASON's video capture system suffers the same problem. Hence this class.
		 */
		try {
		    Thread.sleep(500);  // 1000 milliseconds is one second.
		} catch(InterruptedException ex) {
		    Thread.currentThread().interrupt();
		}

		// Hacked together from mason's Display3D.takeSnapshot() method.
		Simulation sim = (Simulation) state;		
		double time = sim.schedule.getTime();
		if (reenacting)  // See comment above by variable declaration.  
			time = time - Simulation.timeSlice_min;
		String seqNum = String.format("%03d", sequenceNum);
		String fileName = movieDir + "/img_" + seqNum + ".png";

        // Start the image
        display.canvas.beginCapturing(false);
    	// Doesn't seem to actually work
        // It throws an error. java.lang.IllegalStateException: Canvas3D: Not in off-screen mode
        // display.canvas.waitForOffScreenRendering();
		try
        {
			File snapShotFile = new File(fileName);
	        BufferedImage image = display.canvas.getLastImage();
	        String time_fmt = formatTime(time);
	        image = drawTextOnImage(time_fmt, image, 0);
	        PNGEncoder tmpEncoder = new PNGEncoder(image, false,PNGEncoder.FILTER_NONE,9);
	        OutputStream stream = new BufferedOutputStream(new FileOutputStream(snapShotFile));
	        stream.write(tmpEncoder.pngEncode());
	        stream.close();
	        image.flush();  // Just in case -- OS X bug?
        }
		catch (FileNotFoundException e) { } // Fail
		catch (IOException e) { /* Could happen on close? */} // Fail
		sequenceNum ++;
	}

	private BufferedImage drawTextOnImage(String text, BufferedImage image, int space)
	{
	    BufferedImage bi = new BufferedImage(image.getWidth(), image.getHeight() + space, BufferedImage.TRANSLUCENT);
	    Graphics2D g2d = (Graphics2D) bi.createGraphics();
	    g2d.addRenderingHints(new RenderingHints(RenderingHints.KEY_RENDERING, RenderingHints.VALUE_RENDER_QUALITY));
	    g2d.addRenderingHints(new RenderingHints(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON));
	    g2d.addRenderingHints(new RenderingHints(RenderingHints.KEY_TEXT_ANTIALIASING, RenderingHints.VALUE_TEXT_ANTIALIAS_ON));

	    g2d.drawImage(image, 0, 0, null);
	    g2d.setColor(Color.BLACK);
	    g2d.setFont(new Font("Calibri", Font.BOLD, 40));
	    FontMetrics fm = g2d.getFontMetrics();
	    int textWidth = fm.stringWidth(text);

	    //center text at bottom of image in the new space
		int text_x = (int) (0.05 * image.getWidth());
		int text_y = (int) (0.95 * image.getHeight());
	    g2d.drawString(text, text_x, text_y);

	    g2d.dispose();
	    return bi;
	}

	private String formatTime(double minutes)
	{
		int hours = (int) (minutes / 60.); //since both are ints, you get an int
		int min = ((int) minutes) % 60;
		String fmt = String.format("%d:%02d", hours, min);
		return fmt;
	}


//	private void addTime(BufferedImage image)
//	{
//		BufferedImage textImage = new BufferedImage(image.getWidth(), image.getHeight(), image.getType());
//
//		Graphics2D w = (Graphics2D) image.getGraphics();
//		w.drawImage(image, 0, 0, null);
//		int text_x = (int) (0.8 * image.getWidth());
//		int text_y = (int) (0.8 * image.getHeight());
//		String text = "hello";
//		w.drawImage(text, text_x, text_y, null);
//		w.dispose();
//	}




}
