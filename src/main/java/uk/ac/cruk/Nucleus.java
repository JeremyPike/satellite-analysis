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
import ij.plugin.frame.RoiManager;
import net.imagej.ops.OpService;

public class Nucleus {
	
	private SpotCollection satellites;
	private SpotCollection centrosomes;
	private Spot nucleusSpot;
	private double meanDapiIntensity;
	private double meanSatelliteIntensity;
	private double meanCentrosomeIntensity;
	private int numSatellites;
	private int numCentrosomes;
	private double meanNNDistSatCent;
	private double meanNNDistSatSat;
	private boolean isEdge;

	public Nucleus(Spot nucleusSpot, ImagePlus imp, double dapiCentreFrac, double edgeFac) {

		double radius = nucleusSpot.getFeature("RADIUS");
		double radiusPixels = radius / imp.getCalibration().pixelWidth;
		double posX = nucleusSpot.getFeature("POSITION_X");
		double posY = nucleusSpot.getFeature("POSITION_Y");
		OvalRoi oval = new OvalRoi(posX / imp.getCalibration().pixelWidth - radiusPixels / dapiCentreFrac, posY / imp.getCalibration().pixelHeight - radiusPixels / dapiCentreFrac, 2 * radiusPixels / dapiCentreFrac, 2 * radiusPixels / dapiCentreFrac);
		imp.setRoi(oval);
		meanDapiIntensity = imp.getStatistics().mean;
		
		if (posX - 2 * radius <= 0 || posY - 2 * radius <= 0 || posX + 2 * radius >= imp.getWidth() * imp.getCalibration().pixelWidth || posY + 2 * radius >= imp.getHeight() * imp.getCalibration().pixelHeight) {
			isEdge = true;
			nucleusSpot.putFeature("TYPE", 0.0d);
		} else {
			isEdge = false;
			nucleusSpot.putFeature("TYPE", 0.0d);
		}
		
		satellites = new SpotCollection();
		centrosomes = new SpotCollection();
		
		this.nucleusSpot = nucleusSpot;
	}
	
	public SpotCollection getSatellites() {
		return satellites;
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

	public void addSatellite(Spot satellite) {
		satellite.putFeature("TYPE", 1.0d);
		satellites.add(satellite, 0);
	}
	
	public void addCentrosome(Spot centrosome) {
		centrosome.putFeature("TYPE", 2.0d);
		centrosomes.add(centrosome, 0);
	}
	
	public void setNucleusSpotType(double type) {
		nucleusSpot.putFeature("TYPE", type);
	}
	
	public void computeProperties() {
		
		numSatellites = satellites.getNSpots(true);
		numCentrosomes = centrosomes.getNSpots(true);
		
		double[] satelliteIntensities = satellites.collectValues("MEAN_INTENSITY", true);
		meanSatelliteIntensity = 0;
		for (double satInt : satelliteIntensities) {
			meanSatelliteIntensity += satInt;
		}
		meanSatelliteIntensity = meanSatelliteIntensity  / satelliteIntensities.length;
		
		double[] centrosomeIntensities = centrosomes.collectValues("MEAN_INTENSITY", true);
		meanCentrosomeIntensity = 0;
		for (double centInt : centrosomeIntensities) {
			meanCentrosomeIntensity += centInt;
		}
		meanCentrosomeIntensity = meanCentrosomeIntensity  / centrosomeIntensities.length;
		
		Iterator<Spot> satIterator = satellites.iterator(true);
		Spot satellite;
		Spot nearestSat;
		Spot nearestCent;
		meanNNDistSatCent = 0;
		meanNNDistSatSat = 0;
		while (satIterator.hasNext()) {
			satellite = satIterator.next();
			if (satellites.getNSpots(true) > 1) {
				nearestSat = satellites.getNClosestSpots(satellite, 0, 2, true).get(1);
				meanNNDistSatSat += spotSeperation(satellite, nearestSat);
			}
			if (centrosomes.getNSpots(true) != 0) {
				nearestCent = centrosomes.getClosestSpot(satellite, 0, true);
				meanNNDistSatCent += spotSeperation(satellite, nearestCent);
			}
			
		}
		meanNNDistSatCent = meanNNDistSatCent / satellites.getNSpots(true);
		meanNNDistSatSat = meanNNDistSatSat / satellites.getNSpots(true);
	}
	
	public void display(ImagePlus imp) {
		
		SpotCollection allSpots = new SpotCollection();
		allSpots.add(nucleusSpot, 0);
		Iterator<Spot> satIterator = satellites.iterator(true);
		while (satIterator.hasNext()) {
			allSpots.add(satIterator.next(), 0);
		}
		Iterator<Spot> centIterator = centrosomes.iterator(true);
		while (centIterator.hasNext()) {
			allSpots.add(centIterator.next(), 0);
		}
		Model model = new Model();
		model.setSpots(allSpots, false);
		SpotColorGenerator color = new SpotColorGenerator(model);
		color.setFeature("TYPE");
		HyperStackDisplayer displayer =  new HyperStackDisplayer(model, new SelectionModel(model), imp);
		displayer.setDisplaySettings("SPOT_COLORING", color);
		displayer.render();
		displayer.refresh();
	}
	public static double spotSeperation(Spot spot1, Spot spot2) {
		return Math.sqrt(Math.pow(spot1.getFeature("POSITION_X") - spot2.getFeature("POSITION_X"), 2) + Math.pow(spot1.getFeature("POSITION_Y") - spot2.getFeature("POSITION_Y"), 2) + Math.pow(spot1.getFeature("POSITION_Z") - spot2.getFeature("POSITION_Z"), 2));
	}
}
