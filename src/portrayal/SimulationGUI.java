package portrayal;

import java.awt.Color;

import javax.swing.JFrame;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.w3c.dom.Document;
import org.w3c.dom.Node;

import sim.display.Controller;
import sim.display.GUIState;
import sim.display3d.Display3D;
import sim.engine.SimState;
import sim.portrayal3d.continuous.ContinuousPortrayal3D;
import sim.portrayal3d.simple.AxesPortrayal3D;
import sim.portrayal3d.simple.WireFrameBoxPortrayal3D;
import core.Simulation;
import core.SimulationNeutrophil;
import environment.Compartment3D;

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

 * @author Mark N. Read
 *
 */
public abstract class SimulationGUI extends GUIState
{
	// Set to true, this will take pictures of the tissue volume every timestep and save them.
	public static boolean imageVolume = false; 

	// Set to true to prevent items leaving the imaging volume being displayed.
	public static boolean imageOutsideVolume = false;

	// The recording volume can be less than the imaging volume/reachable environment (e.g., thin slice through well).
	// True -> record only a subset of imaging volume. false -> record the entire imaging volume.
	// Currently only works for a slice in the z-axis.
	protected static RecordingVolumeMode imageOutsideRecorded = RecordingVolumeMode.VOLUME;
	protected static enum RecordingVolumeMode {
		VOLUME,  // Portray everything within the imaging volume.
		SLICE,  //  Portray a subset of the imaging volume.
	}

	public static double cellScalar = 1.;  // draw cells bigger, by this factor.

	public Display3D display;
	public JFrame displayFrame;

	public ContinuousPortrayal3D cellPortrayal = new ContinuousPortrayal3D();


	public WireFrameBoxPortrayal3D frame;
	public AxesPortrayal3D axes;

	public SimulationGUI()
	{
		super(new SimulationNeutrophil());
	}

	public SimulationGUI(SimState state)
	{
		super(state);
		try {
			loadParameters(Simulation.parameters);
		} catch(XPathExpressionException e) {
			e.printStackTrace();
		}
	}

	public void start()
	{
		super.start();
		setupPortrayals();
	}

	protected abstract void setupPortrayals();

	public void init(Controller c)
	{
		System.out.println("Motilisim - initialising");
		super.init(c);

		display = new Display3D(1000, 800, this);
		frame = new WireFrameBoxPortrayal3D(
				0.0, 0.0, 0.0,  // One corner, at the origin. Second corner at the opposite extreme.
				Compartment3D.extremeWidth, Compartment3D.extremeHeight, Compartment3D.extremeDepth, 
				Color.BLACK);
		axes = new AxesPortrayal3D(1.0, true);
		display.attach( frame, "frame" );
		display.attach( axes, "axes" );

		// Changes the camera location and focal length.
        display.translate(-Compartment3D.extremeWidth / 2.0,
        				  -Compartment3D.extremeHeight / 2.0,
        				  -Compartment3D.extremeDepth / 2.0);
        display.scale(1.3 / Compartment3D.extremeWidth);
        display.rotateX(180); 	// Flip scene around so camera points at skin.        
        
		displayFrame = display.createFrame();
		displayFrame.setTitle("Tissue Display");
		c.registerFrame(displayFrame);
		displayFrame.setVisible(true);

		display.attach( cellPortrayal, "cells");	
	}

	public void quit()
	{
		super.quit();
		if (displayFrame != null) displayFrame.dispose();
		displayFrame = null;
		display = null;
	}

	public static void loadParameters(Document params) throws XPathExpressionException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath();
		Node n;
		n = (Node) xPath.compile("/params/GUI/cellScalar").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() > 0)  // Tag is present, and not empty.
			cellScalar = Double.parseDouble(n.getTextContent());

		n = (Node) xPath.compile("/params/GUI/visualisationMode").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() > 0)  // Tag is present, and not empty.
		switch(n.getTextContent())
		{
		case "volume": imageOutsideRecorded = RecordingVolumeMode.VOLUME; break;
		case "slice" : imageOutsideRecorded = RecordingVolumeMode.SLICE; break;
		}
		
		n = (Node) xPath.compile("/params/GUI/recordZSlices").evaluate(params, XPathConstants.NODE);
		if (n != null && n.getTextContent().length() > 0)  // Tag is present, and not empty.
			imageVolume = Boolean.parseBoolean(n.getTextContent()); 
		
		// Command line argument can override the parameters file. 
		for(String str: Simulation.cmdLineArgs) {
		    if(str.trim().contains("-m"))  // Take stills to make movie
		    {
		    	imageVolume = true;
		    }
		}
	}
}
