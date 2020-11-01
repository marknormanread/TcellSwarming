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
import org.w3c.dom.Document;
import org.w3c.dom.Node;

import sim.util.Double3D;
import utils.Quaternion;
import core.MigratoryCell;
import core.Simulation;

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
public class Bootstrap implements Orientation, Translation
{
	// source and length of blocks
	private static int blockLength = 0;
	private static String source = "";

	// for managing blocks
	private static HashMap<String, ArrayList<Position>> data = new HashMap<String, ArrayList<Position>>();
	private static ArrayList<BlockStart> blocks = new ArrayList<BlockStart>();
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
	
	public Bootstrap()
	{
		if (fixedSourceTrack)
		{
			ArrayList<String> keyArray = new ArrayList<String>(blocksByTrack.keySet());
			int randomID = Simulation.instance.random.nextInt(keyArray.size());
			sourceTrack = keyArray.get(randomID);
		}
			
		currentBlockIndex = blockLength;  // force new block to be sampled.
		plan();
		// randomize block progress, and re-sample. Avoids artifacts from all cells simultaneously selecting new blocks 
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
		// cell moves along it's x axis. Find its orientation in absolute space, by transforming x-axis. This gives a 
		// unit vector pointing in the direction of the cell's orientation.
		Double3D facing = orientation.transform(MigratoryCell.x_axis);
		// convert unit vector describing cell's orientation in absolute space to a move forward.
		Double3D move = new Double3D(facing.x*length, facing.y*length, facing.z*length);

		return move;
	}

	@Override
	public Quaternion newOrientation(Quaternion orientation, MigratoryCell cell) 
	{
		
		plan();
		// point cell in the same direction as that from which the block is sourced. 
		if (mode == Mode.ORIENTATION)
			orientation = Quaternion.faceVector(currentMovement.displacement);
		else {
			double roll = 0.;
			if (mode == Mode.PITCH)
				roll = Simulation.instance.random.nextDouble() * 2.0 * Math.PI;

			else if (mode == Mode.PITCH_ROLLS)
				roll = this.currentMovement.roll;				
			
			// roll as a quaternion.
			Quaternion rotateQ = Quaternion.representRotation(roll, 
					MigratoryCell.x_axis.x, MigratoryCell.x_axis.y, MigratoryCell.x_axis.z);
			// multiply orientation by rotateQ, because rotateQ is calculated relative to cell, not in absolute space. 
			orientation = orientation.multiply(rotateQ).normalise();  // alter the cell's orientation.

			// Process pitch angle
			final double pitch = this.currentMovement.pitch;
			Quaternion pitchQ = Quaternion.representRotation
					(pitch, MigratoryCell.y_axis.x, MigratoryCell.y_axis.y, MigratoryCell.y_axis.z);
			// multiply orientation by rotateQ, because pitchQ is calculated relative to cell, not in absolute space.
			orientation = orientation.multiply(pitchQ).normalise();
		}
		return orientation;
	}
	
	/** 
	 * Parse source motility data. Data are separated into tracks, and displacement, length and turn angles extracted
	 * for each. Every possible block (given length requirements) in the track is identified and recorded for use by
	 * cells. 
	 * 
	 * @throws IOException
	 */
	private static void initialise() throws IOException 
	{
		List<CSVRecord> records;
		
		// Try to read in the source data file. This can throw exceptions.
		System.out.println(source);
		FileReader reader = new FileReader(source);			
		// Readers headers from file, and skips the header line
		CSVParser parser = CSVFormat.DEFAULT.withHeader().parse(reader);
		records = parser.getRecords();			
				
		// Parse source data, extracting displacements for each track and storing them. 
		CSVRecord record;
		
		// extract positional data 
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
			}
			// Delete any positions that have NaN pitches. 
			Iterator<Position> it = data.get(trackID).iterator();
			while (it.hasNext())
				if (Double.isNaN(it.next().pitch))	
					it.remove();
		}
		
		// Calculate roll angles. 
		if (mode == Mode.PITCH_ROLLS)
		{
			// Add an index (track, starting point) to the possible blocks ArrayList
			for (String trackID : data.keySet())
			{
				ArrayList<Position> positions = data.get(trackID);
				for (int i = 1; i < positions.size() - 1; i++)
				{
					final Double3D u = positions.get(i-1).displacement;
					final Double3D v = positions.get(i).displacement;
					final Double3D w = positions.get(i+1).displacement;
					
					// Following calculations can't be done if there's a zero-length displacement vector
					if (utils.Utils.lengthSq(u) * utils.Utils.lengthSq(v) * utils.Utils.lengthSq(u) != 0.)
					{
						// Projecting U onto V gives u1. u1 + u2 = u, hence we can calculate u2. u2 is the vector lying
			            // on the plane UV that is perpendicular to v.
			            Double3D u1 = utils.Utils.projectVector(u,  v);
			            // Invert this vector, so that it points away from v, not towards it.
			            // Hence, both w2 and u2 will point away from v.
			            Double3D u2 = new Double3D(-1. * (u.x - u1.x),
			            						   -1. * (u.y - u1.y),
			            						   -1. * (u.z - u1.z));
			            // Do the same again for W onto V. 
			            // U2 and w2 will both be perpendicular to v (just as u1 and w1 are parallel to it).
			            Double3D w1 = utils.Utils.projectVector(w, v);			                    
			            Double3D w2 = new Double3D(w.x - w1.x,
			                          			   w.y - w1.y,
			                          			   w.z - w1.z);
			            double angle;
			            if (utils.Utils.lengthSq(u2) * utils.Utils.lengthSq(w2) != 0.0)
			            {
	                     	// Calculate the changes in roll, based on angles between successive orientations.
	                        angle = utils.Utils.rotation_anticlockwise(u2, w2, v, true, true);
			            } else 
			            	/* This can capture some strange boundary conditions. If several displacements have the 
			            	 * exact same orientation, then it's unclear what the roll should be. Consider three 
			            	 * consecutive displacements, if the first two are parallel, and the third is not, it's
			            	 * still unclear with this method what the rotation should be. It will be assigned zero. 
			            	 * Unless all three displacements are non-parallel, zero will be assigned. This shouldn't 
			            	 * occur much in read cell data, but it happens a lot in ballistic motility. */
			            	angle = 0.;
			            positions.get(i).roll = angle;
					}
				}
				
				// Remove any positions for which roll angles haven't been calculated
				Iterator<Position> it = data.get(trackID).iterator();
				while (it.hasNext())
				{
					double roll = it.next().roll;
					if (Double.isNaN(roll))
					{
						it.remove();
					}
				}
			}			
		}
		// Add an index (track, starting point) to the possible blocks ArrayList
		for (String key : data.keySet())
		{			
			ArrayList<Position> positions = data.get(key);
			// Plus one, because a list of 2 displacements and a blockLength of 2 gives 1 usable block 
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
		
		n = (Node) xPath.compile("/params/Motility/Bootstrap/source").evaluate(params, XPathConstants.NODE);
		if (! n.getTextContent().equals(""))  // Don't load if empty. 
			source = n.getTextContent();
		
		initialise();
	}
	
	private static class Position
	{
		Double3D location;
		Double3D displacement;
		double length;
		double pitch = Double.NaN;
		double roll = Double.NaN;
		
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
