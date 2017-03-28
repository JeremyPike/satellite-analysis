package uk.ac.cruk;

import java.util.Iterator;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.Spot;
import fiji.plugin.trackmate.SpotCollection;
import fiji.plugin.trackmate.visualization.SpotColorGenerator;
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer;
import ij.ImagePlus;
import ij.gui.OvalRoi;

/**
 * This class is used to define a single cell and its associated centrosomes and
 * satellites
 * 
 * @author Jeremy Pike
 */
public class Nucleus {

	/**
	 * Collection of {@link Spot} objects referring to satellite detections
	 * associated with this nucleus
	 * 
	 */
	private SpotCollection satellites;
	/**
	 * Collection of satellites after filtering detections close to centrosomes
	 */
	private SpotCollection satellitesFilt;
	/**
	 * Collection of centrosomes associated with is nucleus
	 */
	private SpotCollection centrosomes;
	/**
	 * {@link Spot} for Nucelus Spot
	 */
	private Spot nucleusSpot;
	/**
	 * Mean intensity of Dapi signal measured using a using defined radius about
	 * the centre of nucelusSpot
	 */
	private double meanDapiIntensity;
	/**
	 * Mean intensity across all spots in the filtered satellite collection
	 */
	private double meanSatelliteIntensity;
	/**
	 * Mean intensity across all spots in the centrosome collection
	 */
	private double meanCentrosomeIntensity;
	/**
	 * Number of spots in the filtered satellite collection
	 */
	private int numSatellites;
	/**
	 * Number of centrosomes in the filtered satellite collection
	 */
	private int numCentrosomes;
	/**
	 * Mean distance between each filtered satellite and the nearest centrosome
	 */
	private double meanNNDistSatCent;
	/**
	 * Mean distance between each filtered satellite and the nearest other
	 * filtered satellite
	 */
	private double meanNNDistSatSat;
	/**
	 * Records if centre of nucleusSpot is closer than 2* nuclear radius of the
	 * image edges
	 */
	private boolean isEdge;
	/**
	 * Records if nucleus is mitotic. This decision is based on
	 * meanDapiIntensity
	 * 
	 */
	private boolean isMitotic;

	/**
	 * Construct a new Nucleus and calculates {@link meanDapiIntensity} and
	 * {@link isEdge}
	 * 
	 * @param nucleusSpot
	 *            the {@link Spot} which defines the location and radius of the
	 *            nucleus
	 * 
	 * @param imp
	 *            the image used to calculate the mean dapi intensity
	 * 
	 * @param dapiCentreFrac
	 *            diversion factor for scaling of spot size for mean dapi
	 *            intensity measurement
	 * 
	 * @param edgeFrac
	 *            multiplication factor for nucleus radius to decide if nucleus
	 *            is on the image edge
	 */
	public Nucleus(Spot nucleusSpot, ImagePlus imp, double dapiCentreFrac, double edgeFrac) {

		// set the spot name to empty String
		nucleusSpot.setName(" ");

		// get nucleus radius (should be in microns)
		double radius = nucleusSpot.getFeature("RADIUS");
		// convert to pixels
		double radiusPixels = radius / imp.getCalibration().pixelWidth;
		// get spot centre
		double posX = nucleusSpot.getFeature("POSITION_X");
		double posY = nucleusSpot.getFeature("POSITION_Y");
		// create circular ROI with specified radius at the spot centre
		OvalRoi oval = new OvalRoi(posX / imp.getCalibration().pixelWidth - radiusPixels / dapiCentreFrac,
				posY / imp.getCalibration().pixelHeight - radiusPixels / dapiCentreFrac,
				2 * radiusPixels / dapiCentreFrac, 2 * radiusPixels / dapiCentreFrac);
		// add to dapi image and measure the mean intensity in the ROI
		imp.setRoi(oval);
		meanDapiIntensity = imp.getStatistics().mean;

		// if spot centre is within specified distance of image boundary
		if (posX - edgeFrac * radius <= 0 || posY - edgeFrac * radius <= 0
				|| posX + edgeFrac * radius >= imp.getWidth() * imp.getCalibration().pixelWidth
				|| posY + edgeFrac * radius >= imp.getHeight() * imp.getCalibration().pixelHeight) {
			isEdge = true;
			// set TYPE property, used for spot colouring
			nucleusSpot.putFeature("TYPE", 0.0d);
		} else {
			isEdge = false;
			nucleusSpot.putFeature("TYPE", 2.0d);
		}

		// create empty SpotCollection objects for satellites and centrosomes
		satellites = new SpotCollection();
		centrosomes = new SpotCollection();

		this.nucleusSpot = nucleusSpot;
	}

	/*
	 * Getters and setters
	 */

	public SpotCollection getSatellites() {
		return satellites;
	}

	public SpotCollection getSatellitesFilt() {
		return satellitesFilt;
	}

	public void setSatellites(SpotCollection satellites) {
		this.satellites = satellites;
	}

	public SpotCollection getCentrosomes() {
		return centrosomes;
	}

	public void setCentrosomes(SpotCollection centrosomes) {
		this.centrosomes = centrosomes;
	}

	public boolean getIsEdge() {
		return isEdge;
	}

	public boolean getIsMitotic() {
		return isMitotic;
	}

	public Spot getNucleusSpot() {
		return nucleusSpot;
	}

	public double getMeanDapiIntensity() {
		return meanDapiIntensity;
	}

	public double getMeanSatelliteIntensity() {
		return meanSatelliteIntensity;
	}

	public double getMeanCentrosomeIntensity() {
		return meanCentrosomeIntensity;
	}

	public int getNumSatellites() {
		return numSatellites;
	}

	public int getNumCentrosomes() {
		return numCentrosomes;
	}

	public double getMeanNNDistSatCent() {
		return meanNNDistSatCent;
	}

	public double getMeanNNDistSatSat() {
		return meanNNDistSatSat;
	}

	/**
	 * Adds a {@link Spot} object to the satellite collection
	 * 
	 * @param satellite
	 *            spot to be added
	 */
	public void addSatellite(Spot satellite) {
		// set TYPE for spot colouring
		satellite.putFeature("TYPE", 3.0d);
		// set the spot name to empty String
		satellite.setName(" ");
		satellites.add(satellite, 0);
	}

	/**
	 * Adds a {@link Spot} object to the centrosome collection
	 * 
	 * @param centrosome
	 *            spot to be added
	 */
	public void addCentrosome(Spot centrosome) {
		// set TYPE for spot colouring
		centrosome.putFeature("TYPE", 5.0d);
		// set the spot name to empty String
		centrosome.setName(" ");
		centrosomes.add(centrosome, 0);
	}

	/**
	 * Sets the TYPE for the nucleus spot, used for spot colouring
	 * 
	 * @param type
	 *            TYPE value to be set
	 */
	public void setNucleusSpotType(double type) {
		nucleusSpot.putFeature("TYPE", type);
	}

	/**
	 * Sets the name for the nucleus spot
	 * 
	 * @param name
	 *            name of spot
	 */
	public void setNucleusSpotName(String name) {
		nucleusSpot.setName(name);
	}

	/**
	 * Computes statistics for Fields not set by constructor. This method should
	 * only been called after all relevant satellite and centrosome spots have
	 * been added
	 * 
	 * @param minCentSatDistance
	 *            minimum distance between satellites and centrosomes for
	 *            filtering of satellites
	 * 
	 * @param maxMeanDapiIntensity
	 *            maximum value for mean dapi intensity for nucleus to be
	 *            classed as non-mitotic
	 * 
	 */
	public void computeProperties(double minCentSatDistance, double maxMeanDapiIntensity) {

		Spot satellite;
		Spot nearestCent;
		Spot nearestSat;
		double distSatCent;
		// iterator for all satellites associated with nucleus
		Iterator<Spot> satIterator = satellites.iterator(true);
		// to hold filtered satellites
		satellitesFilt = new SpotCollection();
		while (satIterator.hasNext()) {
			// get each satellite
			satellite = satIterator.next();
			// if there is at least one centrosome associated with this nucleus
			if (centrosomes.getNSpots(true) != 0) {
				// get centrsome which is closest to the satellite
				nearestCent = centrosomes.getClosestSpot(satellite, 0, true);
				// calculate the distance between the nearest centrosome and the
				// satellite
				distSatCent = spotSeperation(satellite, nearestCent, false);
				// if satellite is a sufficient distance from the nearest
				// centrosome add to the filtered collection
				if (distSatCent >= minCentSatDistance) {
					satellitesFilt.add(satellite, 0);
				}
				// if no centrosomes add satellite to the filtered collection
			} else {
				satellitesFilt.add(satellite, 0);
			}
		}

		// get the number of centrosomes and filtered satellites
		numSatellites = satellitesFilt.getNSpots(true);
		numCentrosomes = centrosomes.getNSpots(true);

		// find the mean mean intensity across all filtered satellites
		double[] satelliteIntensities = satellitesFilt.collectValues("MEAN_INTENSITY", true);
		meanSatelliteIntensity = 0;
		for (double satInt : satelliteIntensities) {
			meanSatelliteIntensity += satInt;
		}
		meanSatelliteIntensity = meanSatelliteIntensity / satelliteIntensities.length;

		// find the mean mean intensity across all centrosomes
		double[] centrosomeIntensities = centrosomes.collectValues("MEAN_INTENSITY", true);
		meanCentrosomeIntensity = 0;
		for (double centInt : centrosomeIntensities) {
			meanCentrosomeIntensity += centInt;
		}
		meanCentrosomeIntensity = meanCentrosomeIntensity / centrosomeIntensities.length;

		// iterator for filtered satellites
		satIterator = satellitesFilt.iterator(true);
		// calculate the mean distance across all filtered satellites to the
		// nearest centrosome and filtered satellite
		meanNNDistSatCent = 0;
		meanNNDistSatSat = 0;
		while (satIterator.hasNext()) {
			satellite = satIterator.next();
			if (satellitesFilt.getNSpots(true) > 1) {
				// get closest filtered satellite
				nearestSat = satellitesFilt.getNClosestSpots(satellite, 0, 2, true).get(1);
				// add distance to sum
				meanNNDistSatSat += spotSeperation(satellite, nearestSat, false);
			}
			if (centrosomes.getNSpots(true) != 0) {
				// get closest centrosome
				nearestCent = centrosomes.getClosestSpot(satellite, 0, true);
				// add distance to sum
				meanNNDistSatCent += spotSeperation(satellite, nearestCent, false);
			}

		}
		// divide sum by number of filtered satellites to get mean
		meanNNDistSatCent = meanNNDistSatCent / satellitesFilt.getNSpots(true);
		meanNNDistSatSat = meanNNDistSatSat / satellitesFilt.getNSpots(true);

		// if mean dapi intensity is greater than user specified value then
		// nucleus classified as mitotic
		if (meanDapiIntensity > maxMeanDapiIntensity) {
			isMitotic = true;
			// set TYPE for spot colouring
			nucleusSpot.putFeature("TYPE", 4.0d);
		} else {
			isMitotic = false;
		}

	}

	/**
	 * displays all spots (nucleus centre, satellites and centrsomes) associated
	 * with the Nucelus
	 * 
	 * @param imp
	 *            image used to display detections
	 * 
	 * @param filtered
	 *            if true displays only filtered satellites
	 */
	public void display(ImagePlus imp, boolean filtered) {

		// SpotCollection for all spots associated with nucleus
		SpotCollection allSpots = new SpotCollection();

		// add nucleus spot
		allSpots.add(nucleusSpot, 0);

		// iterate across either all or filtered satellites
		Iterator<Spot> satIterator;
		if (filtered) {
			satIterator = satellitesFilt.iterator(true);
		} else {
			satIterator = satellites.iterator(true);
		}
		while (satIterator.hasNext()) {
			// add to joint collection
			allSpots.add(satIterator.next(), 0);
		}
		// add centrosomes to joint collection
		Iterator<Spot> centIterator = centrosomes.iterator(true);
		while (centIterator.hasNext()) {
			allSpots.add(centIterator.next(), 0);
		}
		// create trackmate model and associate with joint spot collection
		Model model = new Model();
		model.setSpots(allSpots, false);
		// spot colouring based on TYPE property
		SpotColorGenerator color = new SpotColorGenerator(model);
		color.setFeature("TYPE");
		// create and show display
		HyperStackDisplayer displayer = new HyperStackDisplayer(model, new SelectionModel(model), imp);
		displayer.setDisplaySettings(HyperStackDisplayer.KEY_SPOT_COLORING, color);
		displayer.render();
		displayer.refresh();
	}

	/**
	 * calculates the distance between two {@link Spot} objects
	 * 
	 * @param spot1
	 *            the first spot
	 * 
	 * @param spot2
	 *            the second spot
	 * 
	 * @param twoD
	 *            if true only uses the XY coordinates, if false will calculate
	 *            the distance in 3D
	 * 
	 * @return Euclidean distance between the specified spots
	 */
	public static double spotSeperation(Spot spot1, Spot spot2, boolean twoD) {

		// 2D Euclidean distance
		if (twoD) {
			return Math.sqrt(Math.pow(spot1.getFeature("POSITION_X") - spot2.getFeature("POSITION_X"), 2)
					+ Math.pow(spot1.getFeature("POSITION_Y") - spot2.getFeature("POSITION_Y"), 2));
		}
		// 3D Euclidean distance
		else {
			return Math.sqrt(Math.pow(spot1.getFeature("POSITION_X") - spot2.getFeature("POSITION_X"), 2)
					+ Math.pow(spot1.getFeature("POSITION_Y") - spot2.getFeature("POSITION_Y"), 2)
					+ Math.pow(spot1.getFeature("POSITION_Z") - spot2.getFeature("POSITION_Z"), 2));
		}

	}

}
