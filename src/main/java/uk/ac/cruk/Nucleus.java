package uk.ac.cruk;

import java.util.Iterator;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.Spot;
import fiji.plugin.trackmate.SpotCollection;
import fiji.plugin.trackmate.visualization.SpotColorGenerator;
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer;
import ij.ImagePlus;
import net.imagej.ops.OpService;

public class Nucleus {
	
	private SpotCollection satellites;
	private SpotCollection centrosomes;
	private Spot nucleusSpot;
	
	private double meanSatelliteIntensity;
	private double meanCentrosomeIntensity;
	private int numSatellites;
	private int numCentrosomes;
	private double meanNNDistSatCent;
	private double meanNNDistSatSat;
	
	public Nucleus(Spot nucleusSpot) {
		super();
		nucleusSpot.putFeature("TYPE", 0.0d);
		this.nucleusSpot = nucleusSpot;
		satellites = new SpotCollection();
		centrosomes = new SpotCollection();
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
	
	public Spot getNucleusSpot() {
		return nucleusSpot;
	}
	
	public void addSatellite(Spot satellite) {
		satellite.putFeature("TYPE", 1.0d);
		satellites.add(satellite, 0);
	}
	
	public void addCentrosome(Spot centrosome) {
		centrosome.putFeature("TYPE", 2.0d);
		centrosomes.add(centrosome, 0);
	}
	
	
	public void computeStats() {
		
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
	
}
