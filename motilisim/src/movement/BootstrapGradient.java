package movement;

import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.commons.math3.exception.MathArithmeticException;
import org.apache.commons.math3.geometry.euclidean.threed.Rotation;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.w3c.dom.Document;
import org.w3c.dom.Node;

import core.MigratoryCell;
import core.Simulation;
import environment.BoundedCylinder;
import sim.util.Double3D;
import utils.Quaternion;

/**
 * Rather than modelling motility, this Class re-samples an existing data set's motility, recreating it here in 
 * simulation. A source file must be supplied. The profile of turns and jumps (speeds) are then extracted. This is done
 * in blocks (length specified by user). The blocks allow temporal patterns and autocorrelations to be maintained in
 * the recreated data. 
 * 
 * There is no preservation of heterogeneity across tracks, beyond the block sampling. 
 * 
 * 
 * @author Mark N. Read, 2018
 *
 */
public class BootstrapGradient implements Orientation, Translation
{
	// source and length of blocks
	private static int blockLength = 0;
	private static String source = "";
	private static Vector3D source_gradient;
	
	// for managing blocks
	private static HashMap<String, ArrayList<Position>> data = new HashMap<String, ArrayList<Position>>();
	// Used as index into data if agents are permitted to sample blocks from any tracks. 
	private static ArrayList<BlockStart> blocks = new ArrayList<BlockStart>();
	// Used as index into data if agents are restricted to bootstrap from a single track (and maintain this over time)
	private static HashMap<String, ArrayList<BlockStart>> blocksByTrack = new HashMap<String, ArrayList<BlockStart>>();
	
	// pertaining to individual cells
	private BlockStart currentBlock;
	private int currentBlockIndex = 0; // each cell must retain knowledge of where it is in the selected block.
	private Position currentMovement = null;
	
	private static enum Mode
	{
		PITCH,
		PITCH_ROLLS,
		ORIENTATION
	};
	private static Mode mode;
	// if true, simulated walkers stick to a single source cell for selecting blocks
	private static boolean fixedSourceTrack = false;
	private String sourceTrack;
	
	final static Vector3D x = new Vector3D(1, 0, 0);
	final static Vector3D y = new Vector3D(0, 1, 0);
	final static Vector3D z = new Vector3D(0, 0, 1);
	
	public BootstrapGradient()
	{
		if (fixedSourceTrack)
		{
			ArrayList<String> keyArray = new ArrayList<String>(blocksByTrack.keySet());
			int randomID = Simulation.instance.random.nextInt(keyArray.size());
			sourceTrack = keyArray.get(randomID);
		}
			
		currentBlockIndex = blockLength;  // force new block to be sampled.
		plan();
		// Randomize block progress, and re-sample. Avoids artifacts from all cells simultaneously selecting new blocks 
		currentBlockIndex = Simulation.instance.random.nextInt(blockLength);
		plan();
		
	}
	
	/**  Sample the next real-step  
	 */
	private void plan()
	{		
		if (currentBlockIndex >= blockLength)  // reached end of block, resample.
		{
			if (fixedSourceTrack)
			{
				final int newBlockIndex = Simulation.instance.random.nextInt(blocksByTrack.get(sourceTrack).size());
				currentBlock = blocksByTrack.get(sourceTrack).get(newBlockIndex);
			} else {
				// Blocks can be sourced form any track. 
				final int newBlockIndex = Simulation.instance.random.nextInt(blocks.size());			
				currentBlock = blocks.get(newBlockIndex);
			}
			currentBlockIndex = 0;
		}
		currentMovement = data.get(currentBlock.trackID).get(currentBlock.start + currentBlockIndex);		
		currentBlockIndex++;
	}
	
	@Override
	public Double3D move(Quaternion orientation) 
	{		
		final double length = this.currentMovement.length;	
		// Cell moves along it's x axis. 
		// Find its orientation in absolute space, by transforming cell's x-axis.  
		// This gives a unit vector pointing in the direction of the cell's orientation in absolute space.		
		// Translate movement along cell's axis into a displacement in absolute space. 
		final Double3D move = orientation.transform(MigratoryCell.x_axis.multiply(length));
			
		return move;
	}

	@Override
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell) 
	{
		plan();  // Maintenance on block sampling. 
		// Direction of gradient at cell's current location, relative to the cell.
		Double3D grad = cell.perceiveAttractant().gradientDirection();
		grad = orientation.transform(grad);  // Transform to a gradient direction in absolute space
				
		Vector3D presentGradient = new Vector3D(grad.x, grad.y, grad.z);		
		Rotation r__ = null;
		
		if (presentGradient.dotProduct(z) == 1) // Vectors are parallel; almost never happens.
		{
			// If presentGradient parallel to y, no notion of "up can be formed. There is no g__y plane
			r__ = new Rotation(x, presentGradient);  
		} else {
			// If angle between presentGradient and y is not same as angle between x and y, then the half plane lying
			// on presentGradientY is used instead. 
			try {
				r__ = new Rotation(x, z, presentGradient, z);
			} catch(MathArithmeticException e) {
				// Precision error near gradient vector length = 0. Assign random direction to induce undirected motion.
				presentGradient = new Vector3D(Simulation.instance.random.nextDouble()-0.5,
											   Simulation.instance.random.nextDouble()-0.5,
											   Simulation.instance.random.nextDouble()-0.5);
				r__ = new Rotation(x, z, presentGradient, z);				
			} catch(Exception e) { 
				System.out.println("Fatal error in chemokine gradient based bootstrapping. ");
				e.printStackTrace();
			}
		}
		
		Vector3D d__ = r__.applyTo(currentMovement.displacementWRTGradient);
		Double3D displacement = new Double3D(d__.getX(), d__.getY(), d__.getZ());
		orientation = Quaternion.faceVector(displacement);
		
		return orientation;
	}
	
	/** Parse source motility data. Data are separated into tracks, and displacement, length and turn angles extracted
	 * for each. Every possible block (given length requirements) in the track is identified and recorded for use by
	 * cells. 
	 * 
	 * @throws IOException
	 */
	private static void initialise() throws IOException 
	{
		List<CSVRecord> records;
		
		// try to read in the source data file. This can throw exceptions.
		System.out.println(source);
		FileReader reader = new FileReader(source);			
		// readers headers from file, and skips the header line
		CSVParser parser = CSVFormat.DEFAULT.withHeader().parse(reader);
		records = parser.getRecords();			
				
		// parse source data, extracting displacements for each track and storing them. 
		CSVRecord record;
		
		// extract positional data, one line in the source file at a time. 
		for (int i = 0; i < records.size(); i++)  
		{
			record = records.get(i);
			String recID = record.get("TrackID");
			Double3D currLoc = new Double3D(Double.parseDouble(record.get("Position_X")),
										    Double.parseDouble(record.get("Position_Y")),
										    Double.parseDouble(record.get("Position_Z")));
			
			if (! data.containsKey(recID))  // initialize if first time track ID is encountered
				data.put(recID, new ArrayList<Position>());
				
			data.get(recID).add(new Position(currLoc));
		}
		
		// Calculate pitch angles, displacements vectors and speeds
		for (String trackID : data.keySet())
		{
			ArrayList<Position> positions = data.get(trackID);
			for (int i = 1; i < positions.size() - 1; i++)
			{
				final Double3D a = positions.get(i-1).location;
				final Double3D b = positions.get(i).location;
				final Double3D c = positions.get(i+1).location;
				
				final Double3D disp1 = b.subtract(a);
				final Double3D disp2 = c.subtract(b);
				
				positions.get(i).displacement = disp1;
				final double length1 = utils.Utils.length(disp1);
				positions.get(i).length = length1;
				
				if (length1 * utils.Utils.lengthSq(disp2) != 0.)
				{
					final double turn = utils.Utils.angleBetweenVectors(disp1, disp2);
					positions.get(i).pitch = turn;					
				}
				
				/* calculate displacement w.r.t. the chemoattractant gradient
				 * The rotations defined here attempt to preserve orientation of by maintaining a notion of "up" (though any
				 * direction will do, we use the y axis). 
				 * This is needed to preserve roll autocorrelation. 
				 * Without the "up" preservation, the rotation from x axis to gradient can reorientatate y and z - there are 
				 * infinite rotations between two vectors in 3D. 
				 * The resultant reorientations of x and source gradient may not be the same as the subsequent rotation between
				 * x and the current gradient. 
				 * If this discrepancy keeps occurring, autocorrelation will be destroyed. This doesn't matter in the extreme
				 * case where a source cell follows the source gradient exactly. 
				 * There's a special case for when the gradient aligns perfectly with the vector being considered "up", as 
				 * the rotation becomes under-specified. */
				Rotation r;
				if (source_gradient.dotProduct(y) == 1)
					// vectors are parallel
					r = new Rotation(x, source_gradient);
				else
					r = new Rotation(x, z, source_gradient, z);	
				
				// displacement with respect to the gradient
				// conversions between Double3D and Vector3D, Rotations (apache commons math) uses Vector3D
				Vector3D d_ = r.applyInverseTo(new Vector3D(disp1.x, disp1.y, disp1.z));
				
				positions.get(i).displacementWRTGradient = d_;
			}
			// delete any positions that have NaN pitches. 
			Iterator<Position> it = data.get(trackID).iterator();
			while (it.hasNext())
				if (it.next().displacementWRTGradient == null)	
					it.remove();
		}
		
		// add an index (track, starting point) to the possible blocks ArrayList
		for (String key : data.keySet())
		{			
			ArrayList<Position> positions = data.get(key);
			// plus one, because a list of 2 displacements and a blockLength of 2 gives 1 usable block 
			final int extractableBlocks = positions.size() - blockLength + 1;
			for (int i = 0; i < extractableBlocks; i++)
			{
				if (! blocksByTrack.containsKey(key))  // If this track isn't already registered, initialise
					blocksByTrack.put(key, new ArrayList<BlockStart>());
				
				BlockStart bs = new BlockStart(key, i); 
				blocksByTrack.get(key).add(bs);
				blocks.add(bs);
			}
		}
	}
	
	public static void loadParameters(Document params) throws XPathExpressionException, IOException
	{
		XPath xPath =  XPathFactory.newInstance().newXPath(); 
		Node n;
		
		n = (Node) xPath.compile("/params/Motility/Bootstrap/blockLength").evaluate(params, XPathConstants.NODE);
		blockLength = Integer.parseInt(n.getTextContent());

		// select mode of bootstrapping
		n = (Node) xPath.compile("/params/Motility/Bootstrap/mode").evaluate(params, XPathConstants.NODE);
		if (n.getTextContent().equals("orientation"))
			mode = Mode.ORIENTATION;
		else if (n.getTextContent().equals("pitch"))
			mode = Mode.PITCH;
		else if (n.getTextContent().equals("pitchRoll"))
			mode = Mode.PITCH_ROLLS;
		
		n = (Node) xPath.compile("/params/Motility/Bootstrap/fixedSourceTrack").evaluate(params, XPathConstants.NODE);
		fixedSourceTrack = Boolean.parseBoolean(n.getTextContent());
		
		n = (Node) xPath.compile("/params/Motility/BootstrapGradient/source").evaluate(params, XPathConstants.NODE);
		if (! n.getTextContent().equals(""))  // don't load if empty. 
			source = n.getTextContent();
		n = (Node) xPath.compile("/params/Motility/BootstrapGradient/gradient_x").evaluate(params, XPathConstants.NODE);
		double gradient_x = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/BootstrapGradient/gradient_y").evaluate(params, XPathConstants.NODE);
		double gradient_y = Double.parseDouble(n.getTextContent());
		n = (Node) xPath.compile("/params/Motility/BootstrapGradient/gradient_z").evaluate(params, XPathConstants.NODE);
		double gradient_z = Double.parseDouble(n.getTextContent());
		source_gradient = new Vector3D(gradient_x, gradient_y, gradient_z);
		initialise();
	}
	
	private static class Position
	{
		Double3D location;
		Double3D displacement = null;
		/* Displacement with respect to the gradient. This displacement is expressed through basis vectors i, j, k, 
		 * which is a transformation of absolute spatial axes x, y, and z. Component i is parallel to the gradient 
		 * direction at the cell's current location. Component j lies in the plane iy, hence it always points towards y.
		 * This is needed to ensure that ijk bases maintain a relatively constant orientation; "up" is always "up", and
		 * this helps ensure that roll autocorrelations are maintained. k lies perpendicular to the plane ij. */
		Vector3D displacementWRTGradient = null; 
		double length;
		double pitch = Double.NaN;
		
		
		public Position(Double3D location)
		{	
			this.location = location;
		}
	}
	
	private static class BlockStart
	{
		String trackID;
		int start;
		
		public BlockStart(String tid, int s)
		{
			this.trackID = tid;
			this.start = s;
		}
	}
}

	