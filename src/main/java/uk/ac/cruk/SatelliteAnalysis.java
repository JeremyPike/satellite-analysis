
package uk.ac.cruk;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.Future;

import org.scijava.command.Command;
import org.scijava.command.CommandModule;
import org.scijava.module.ModuleItem;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.Spot;
import fiji.plugin.trackmate.SpotCollection;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.detection.DownsampleLogDetectorFactory;
import fiji.plugin.trackmate.detection.LogDetectorFactory;
import fiji.plugin.trackmate.features.FeatureFilter;
import fiji.plugin.trackmate.features.spot.SpotIntensityAnalyzerFactory;
import fiji.plugin.trackmate.visualization.SpotColorGenerator;
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer;
import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.ResultsTable;
import ij.process.LUT;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.legacy.LegacyService;
import net.imagej.ops.Op;
import net.imagej.ops.OpService;
import net.imglib2.FinalInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

/**
 * This is a imageJ command for spatial analysis of nuclei, centrosomes and
 * satellites.
 * 
 */

@Plugin(type = Command.class, menuPath = "Plugins>Satellite Analysis")
public class SatelliteAnalysis<T extends RealType<T>> implements Command {

	@Parameter
	private Dataset currentData;

	@Parameter
	private ImageJ ij;

	@Parameter
	private UIService uiService;

	@Parameter
	private OpService opService;

	@Parameter
	private DatasetService datasetService;

	@Parameter(label = "Radius for Gaussian blur (DAPI channel) (microns): ", persist = false, min = "0", max = "10")
	private double dapiSigmaMicrons = 1;

	@Parameter(label = "Radius for nuclei spot detection (microns): ", persist = false, min = "1", max = "20")
	private double dapiSpotRadius = 7;

	@Parameter(label = "Quality threshold for nuclei spots: ", persist = false, min = "0")
	private double dapiSpotQualityThresh = 0.075;

	@Parameter(label = "Radius for satellite spot detection (microns): ", persist = false, min = "0.01", max = "20")
	private double satSpotRadius = 0.2;

	@Parameter(label = "Mean intenisty threshold for satelllite spots: ", persist = false, min = "0")
	private double satSpotMeanIntThresh = 26;

	@Parameter(label = "Radius for Gaussian blur (Centrosome channel) (microns): ", persist = false, min = "0", max = "10")
	private double centSigmaMicrons = 0.2;

	@Parameter(label = "Radius for rolling ball BS (Centrosome channel) (microns): ", persist = false, min = "0", max = "10")
	private double centRollingBallRad = 2;

	@Parameter(label = "Radius for centrosome spot detection (microns): ", persist = false, min = "0.01", max = "20")
	private double centSpotRadius = 0.75;

	@Parameter(label = "Max intenisty threshold for centrosome spots: ", persist = false, min = "0")
	private double centSpotMaxIntThresh = 18;

	@Parameter(label = "Division factor for measurement within DAPI spots ", persist = false, min = "1")
	private double dapiCentreFrac = 4;

	@Parameter(label = "Multiplication edge factor for minimum distance from image edge for Dapi spots ", persist = false, min = "0.1")
	private double edgeFac = 2;

	@Parameter(label = "Multiplication factor for median to calculate maximum nuceli medan intensity ", persist = false, min = "0.1")
	private double nuceliMedianFac = 1.5;

	@Parameter(label = "Filter satellites closer to centrosomes than centrosome spot radius", persist = false)
	private boolean filterSatellites = true;

	@Parameter(label = "Display spot detection overlays on maximal projection", persist = false)
	private boolean displayMPOverlay = true;

	@Override
	public void run() {
		// Create legacy service
		LegacyService legacy = ij.get(LegacyService.class);

		// Get series name
		String seriesName = currentData.getName();

		// Get IJ2 Img object
		final Img<T> image = (Img<T>) currentData.getImgPlus();

		// Pull out data size
		long width = image.dimension(0);
		long height = image.dimension(1);
		long numChannels = image.dimension(2);
		long numSlices = image.dimension(3);

		// Get Dapi data
		FinalInterval slice = FinalInterval.createMinSize(0, 0, 0, 0, width, height, 1, numSlices);
		RandomAccessibleInterval<T> dapiImage = opService.transform().crop(image, slice);

		// Dimensions for max proj image
		long[] projectedDimensions = new long[] { width, height };
		// Generate memory for max proj with type same as input data
		Img<T> dapiImageMP = (Img<T>) opService.run("create.img", projectedDimensions, image.firstElement());
		// Maximal projection along z axis
		Op maxProj = opService.op("stats.max", dapiImage);
		opService.run("transform.project", dapiImageMP, dapiImage, opService.op("stats.max", dapiImage), 2);

		// Convert image to FloatType for better numeric precision
		Img<FloatType> dapiImageMPBlur = opService.convert().float32(dapiImageMP);
		// Gaussian filtering
		opService.filter().gauss(dapiImageMPBlur, dapiImageMPBlur, dapiSigmaMicrons / currentData.averageScale(0));

		// Convert to IJ1 ImagePlus for Trackmate
		ImagePlus dapiImageMPBlurImp = legacy.getImageMap().registerDataset(datasetService.create(dapiImageMPBlur));
		ImagePlus dapiImageMPImp = legacy.getImageMap().registerDataset(datasetService.create(dapiImageMP));
		IJ.run(dapiImageMPBlurImp, "Properties...", "channels=1 slices=1 frames=1 unit=micron pixel_width="
				+ currentData.averageScale(0) + " pixel_height=" + currentData.averageScale(1) + " voxel_depth=1");
		IJ.run(dapiImageMPImp, "Properties...", "channels=1 slices=1 frames=1 unit=micron pixel_width="
				+ currentData.averageScale(0) + " pixel_height=" + currentData.averageScale(1) + " voxel_depth=1");

		// Use Trackmate to count the number of nuclei
		Settings dapiDettings = new Settings();
		dapiDettings.detectorFactory = new DownsampleLogDetectorFactory();
		dapiDettings.addSpotFilter(new FeatureFilter("QUALITY", dapiSpotQualityThresh, true));
		Map<String, Object> map = dapiDettings.detectorFactory.getDefaultSettings();
		map.put("RADIUS", dapiSpotRadius);
		map.put("DO_MEDIAN_FILTERING", false);
		map.put("DOWNSAMPLE_FACTOR", 4);
		map.put("DO_SUBPIXEL_LOCALIZATION", false);
		map.put("TARGET_CHANNEL", 1.0d);
		map.put("THRESHOLD", 0.0d);
		dapiDettings.detectorSettings = map;
		dapiDettings.setFrom(dapiImageMPBlurImp);

		SpotCollection dapiSpots = countSpotsTrackmate(dapiDettings, false);

		// Get satellite data
		slice = FinalInterval.createMinSize(0, 0, 1, 0, width, height, 1, numSlices);
		RandomAccessibleInterval<T> satImage = opService.transform().crop(image, slice);
		// Convert image to FloatType for better numeric precision
		Img<FloatType> satImageFloat = (Img<FloatType>) opService.run("convert.float32", satImage);

		// Convert to IJ1 ImagePlus for Trackmate
		ImagePlus satImageImp = legacy.getImageMap().registerDataset(datasetService.create(satImageFloat));
		IJ.run(satImageImp, "Properties...",
				"channels=1 slices=" + image.dimension(3) + " frames=" + image.dimension(4)
						+ " unit=micron pixel_width=" + currentData.averageScale(0) + " pixel_height="
						+ currentData.averageScale(1) + " voxel_depth=" + currentData.averageScale(3));

		// Use Trackmate to count the number of nuclei
		Settings satSettings = new Settings();
		satSettings.detectorFactory = new LogDetectorFactory();
		satSettings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<>());
		satSettings.addSpotFilter(new FeatureFilter("MEAN_INTENSITY", satSpotMeanIntThresh, true));
		map = satSettings.detectorFactory.getDefaultSettings();
		map.put("RADIUS", satSpotRadius);
		map.put("DO_MEDIAN_FILTERING", false);
		map.put("DO_SUBPIXEL_LOCALIZATION", true);
		map.put("TARGET_CHANNEL", 1.0d);
		map.put("THRESHOLD", 0.0d);
		satSettings.detectorSettings = map;
		satSettings.setFrom(satImageImp);

		SpotCollection satSpots = countSpotsTrackmate(satSettings, false);

		// Get centrosome data
		slice = FinalInterval.createMinSize(0, 0, 2, 0, width, height, 1, numSlices);
		RandomAccessibleInterval<T> centImage = opService.transform().crop(image, slice);
		// Convert image to FloatType for better numeric precision
		Img<FloatType> centImageBlur = (Img<FloatType>) opService.run("convert.float32", centImage);

		// Img<FloatType> centImageBlur =
		// opService.convert().float32((IterableInterval<C>) centImage);

		// Gaussian filtering
		opService.filter().gauss(centImageBlur, centImageBlur, centSigmaMicrons / currentData.averageScale(0),
				centSigmaMicrons / currentData.averageScale(1), centSigmaMicrons / currentData.averageScale(3));

		// Convert to IJ1 ImagePlus for Trackmate
		ImagePlus centImageBlurImp = legacy.getImageMap().registerDataset(datasetService.create(centImageBlur));
		ImagePlus centImageImp = legacy.getImageMap().registerDataset(datasetService.create(centImage));

		IJ.run(centImageBlurImp, "Properties...",
				"channels=1 slices=" + image.dimension(3) + " frames=" + image.dimension(4)
						+ " unit=micron pixel_width=" + currentData.averageScale(0) + " pixel_height="
						+ currentData.averageScale(1) + " voxel_depth=" + currentData.averageScale(3));
		IJ.run(centImageImp, "Properties...",
				"channels=1 slices=" + image.dimension(3) + " frames=" + image.dimension(4)
						+ " unit=micron pixel_width=" + currentData.averageScale(0) + " pixel_height="
						+ currentData.averageScale(1) + " voxel_depth=" + currentData.averageScale(3));
		// Rolling ball background subtraction
		IJ.run(centImageBlurImp, "Subtract Background...",
				"rolling=" + centRollingBallRad / currentData.averageScale(0) + " stack");

		// Use Trackmate to count the number of nuclei
		Settings centSettings = new Settings();
		centSettings.detectorFactory = new LogDetectorFactory();
		centSettings.addSpotAnalyzerFactory(new SpotIntensityAnalyzerFactory<>());
		centSettings.addSpotFilter(new FeatureFilter("MAX_INTENSITY", centSpotMaxIntThresh, true));
		map = centSettings.detectorFactory.getDefaultSettings();
		map.put("RADIUS", centSpotRadius);
		map.put("DO_MEDIAN_FILTERING", false);
		map.put("DO_SUBPIXEL_LOCALIZATION", true);
		map.put("TARGET_CHANNEL", 1.0d);
		map.put("THRESHOLD", 0.0d);
		centSettings.detectorSettings = map;
		centSettings.setFrom(centImageBlurImp);

		SpotCollection centSpots = countSpotsTrackmate(centSettings, centImageImp, false);

		int numUnfilteredNuclei = dapiSpots.getNSpots(true);
		Iterator<Spot> dapiIterator = dapiSpots.iterator(true);
		Nucleus[] nucleiUnfiltered = new Nucleus[numUnfilteredNuclei];
		double[] nuceliMean = new double[numUnfilteredNuclei];

		for (int i = 0; i < numUnfilteredNuclei; i++) {
			Spot dapiSpot = dapiIterator.next();
			nucleiUnfiltered[i] = new Nucleus(dapiSpot, dapiImageMPImp, dapiCentreFrac, edgeFac);
			nuceliMean[i] = nucleiUnfiltered[i].getMeanDapiIntensity();
		}

		double nucleiMeanMedian;
		Arrays.sort(nuceliMean);
		if (nuceliMean.length % 2 == 0) {
			nucleiMeanMedian = (nuceliMean[nuceliMean.length / 2] + nuceliMean[nuceliMean.length / 2 - 1]) / 2;
		} else {
			nucleiMeanMedian = nuceliMean[nuceliMean.length / 2];
		}

		ArrayList<Nucleus> nucleiAL = new ArrayList<Nucleus>();

		for (int i = 0; i < numUnfilteredNuclei; i++) {
			boolean nucCheck = true;
			for (int j = 0; j < numUnfilteredNuclei; j++) {
				double nucSepDist = Nucleus.spotSeperation(nucleiUnfiltered[i].getNucleusSpot(),
						nucleiUnfiltered[j].getNucleusSpot(), true);
				if (i != j && nucSepDist <= dapiSpotRadius
						&& nucleiUnfiltered[i].getMeanDapiIntensity() < nucleiUnfiltered[j].getMeanDapiIntensity()) {
					nucCheck = false;

					nucleiUnfiltered[i].setNucleusSpotType(4.0d);
				}

			}
			if (nucCheck) {
				nucleiAL.add(nucleiUnfiltered[i]);
			}
		}

		Nucleus[] nuclei = nucleiAL.toArray(new Nucleus[nucleiAL.size()]);
		int numNuclei = nuclei.length;

		Iterator<Spot> centIterator = centSpots.iterator(true);
		Spot cent;
		double minDist;
		int minInd;
		while (centIterator.hasNext()) {
			cent = centIterator.next();
			minDist = Nucleus.spotSeperation(cent, nuclei[0].getNucleusSpot(), true);
			minInd = 0;
			for (int i = 1; i < numNuclei; i++) {
				double dist = Nucleus.spotSeperation(cent, nuclei[i].getNucleusSpot(), true);
				if (dist < minDist) {
					minDist = dist;
					minInd = i;
				}
			}
			nuclei[minInd].addCentrosome(cent);
		}

		Iterator<Spot> satIterator = satSpots.iterator(true);
		Spot sat;

		while (satIterator.hasNext()) {
			sat = satIterator.next();
			minDist = Nucleus.spotSeperation(sat, nuclei[0].getNucleusSpot(), true);
			minInd = 0;
			for (int i = 1; i < numNuclei; i++) {
				double dist = Nucleus.spotSeperation(sat, nuclei[i].getNucleusSpot(), true);
				if (dist < minDist) {
					minDist = dist;
					minInd = i;
				}
			}
			nuclei[minInd].addSatellite(sat);
		}

		ResultsTable rt = ResultsTable.getResultsTable();

		for (int i = 0; i < numNuclei; i++) {

			if (filterSatellites) {
				nuclei[i].computeProperties(centSpotRadius, nucleiMeanMedian * nuceliMedianFac);
			} else {
				nuclei[i].computeProperties(0.0d, nucleiMeanMedian * nuceliMedianFac);
			}
			if (!nuclei[i].getIsEdge() && !nuclei[i].getIsMitotic()) {
				rt.incrementCounter();
				rt.addValue("Series name", seriesName);
				rt.addValue("Mean DAPI intensity", nuclei[i].getMeanDapiIntensity());
				rt.addValue("Number of satellites", nuclei[i].getNumSatellites());
				rt.addValue("Number of centrosomes", nuclei[i].getNumCentrosomes());
				rt.addValue("Mean satellite intensity", nuclei[i].getMeanSatelliteIntensity());
				rt.addValue("Mean centrosome intensity", nuclei[i].getMeanCentrosomeIntensity());
				rt.addValue("Mean NN distance satellite-satellite", nuclei[i].getMeanNNDistSatSat());
				rt.addValue("Mean NN distance satellite-centrosome", nuclei[i].getMeanNNDistSatCent());
			}

		}
		rt.show("Results");

		if (displayMPOverlay) {
			// Dimensions for max proj image
			projectedDimensions = new long[] { width, height, numChannels };
			// Generate memory for max proj with type same as input data
			Img<T> allImageMP = (Img<T>) opService.run("create.img", projectedDimensions, image.firstElement());
			// Maximal projection along z axis
			maxProj = opService.op("stats.max", currentData);
			opService.run("transform.project", allImageMP, currentData, opService.op("stats.max", currentData), 3);
			ImagePlus allImageMPImp = legacy.getImageMap().registerDataset(datasetService.create(allImageMP));
			IJ.run(allImageMPImp, "Properties...",
					"channels=" + numChannels + " slices=1 frames=1 unit=micron pixel_width="
							+ currentData.averageScale(0) + " pixel_height=" + currentData.averageScale(1)
							+ " voxel_depth=1");
			CompositeImage impComp = new CompositeImage(allImageMPImp);
			impComp.setChannelLut(LUT.createLutFromColor(Color.BLUE), 1);
			impComp.setChannelLut(LUT.createLutFromColor(Color.GREEN), 2);
			impComp.setChannelLut(LUT.createLutFromColor(Color.RED), 3);

			allImageMPImp = (ImagePlus) impComp.clone();
			allImageMPImp.setDisplayMode(IJ.COMPOSITE);

			displayAll(allImageMPImp, nuclei, true);

		}

	}

	/**
	 * This main function serves for development purposes.
	 *
	 * @param args
	 *            whatever, it's ignored
	 * @throws Exception
	 */
	public static void main(final String... args) throws Exception {

		final ImageJ ij = new ImageJ();
		ij.ui().showUI();
		String id = "C:\\Users\\pike01\\Documents\\testSatSeries.tif";
		Dataset dataset = (Dataset) ij.io().open(id);
		ij.ui().show(dataset);
		//ij.command().run(SatelliteAnalysis.class, true);
	
		Map<String, Object> pluginSettings = new HashMap<String, Object>();
		Iterator<ModuleItem<?>> inputs = ij.command().getCommand(SatelliteAnalysis.class).inputs().iterator();
		while (inputs.hasNext()) {
			ModuleItem<?> input = inputs.next();
			if (input.getDefaultValue() != null) 
			pluginSettings.put(input.getName(), input.getDefaultValue());
		}
		
	
	}

	private static SpotCollection countSpotsTrackmate(Settings settings, boolean display) {

		Model model = new fiji.plugin.trackmate.Model();
		TrackMate trackmate = new TrackMate(model, settings);
		// Check input is ok
		boolean ok = trackmate.checkInput();
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}
		// Find spots
		ok = trackmate.execDetection();
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}
		// Compute spot features
		ok = trackmate.computeSpotFeatures(true);
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}

		// Filter spots
		ok = trackmate.execSpotFiltering(true);
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}

		if (display) {
			SelectionModel selectionModel = new SelectionModel(model);
			HyperStackDisplayer displayer = new HyperStackDisplayer(model, selectionModel, settings.imp);
			displayer.render();
			displayer.refresh();
		}
		// Return spot collection
		return model.getSpots();
	}

	private static SpotCollection countSpotsTrackmate(Settings settings, ImagePlus imp, boolean display) {

		Model model = new fiji.plugin.trackmate.Model();
		TrackMate trackmate = new TrackMate(model, settings);
		// Check input is ok
		boolean ok = trackmate.checkInput();
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}
		// Find spots
		ok = trackmate.execDetection();
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}
		// Compute spot features
		ok = trackmate.computeSpotFeatures(true);
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}

		// Filter spots
		ok = trackmate.execSpotFiltering(true);
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}

		if (display) {
			SelectionModel selectionModel = new SelectionModel(model);
			HyperStackDisplayer displayer = new HyperStackDisplayer(model, selectionModel, settings.imp);
			displayer.render();
			displayer.refresh();
		}

		// Compute spot features with imp
		settings.setFrom(imp);
		ok = trackmate.computeSpotFeatures(true);
		if (ok == false) {
			System.out.println(trackmate.getErrorMessage());
		}
		// Return spot collection
		return model.getSpots();
	}

	private void displayAll(ImagePlus imp, Nucleus[] nuclei, boolean filtered) {

		SpotCollection allSpots = new SpotCollection();
		Spot spot;
		SpotCollection satellites;
		SpotCollection centrosomes;
		Iterator<Spot> satIterator;
		Iterator<Spot> centIterator;
		for (int i = 0; i < nuclei.length; i++) {

			if (!filtered || !nuclei[i].getIsEdge() && !nuclei[i].getIsMitotic()) {
				spot = nuclei[i].getNucleusSpot();
				spot.putFeature("POSITION_Z", 0.0d);
				allSpots.add(spot, 0);
				if (filtered) {
					satellites = nuclei[i].getSatellitesFilt();
				} else {
					satellites = nuclei[i].getSatellites();
				}
				centrosomes = nuclei[i].getCentrosomes();
				satIterator = satellites.iterator(true);
				while (satIterator.hasNext()) {
					spot = satIterator.next();
					spot.putFeature("POSITION_Z", 0.0d);
					allSpots.add(spot, 0);
				}
				centIterator = centrosomes.iterator(true);
				while (centIterator.hasNext()) {
					spot = centIterator.next();
					spot.putFeature("POSITION_Z", 0.0d);
					allSpots.add(spot, 0);
				}

			}

		}

		Model modelAll = new Model();
		modelAll.setSpots(allSpots, false);
		SpotColorGenerator color = new SpotColorGenerator(modelAll);
		color.setFeature("TYPE");

		Iterator<Spot> allIterator = allSpots.iterator(true);

		while (allIterator.hasNext()) {
			Spot tSpot = allIterator.next();
			System.out.println(tSpot.getFeature("TYPE"));
			System.out.println(color.color(tSpot).toString());
		}

		HyperStackDisplayer displayer = new HyperStackDisplayer(modelAll, new SelectionModel(modelAll), imp);
		
		displayer.setDisplaySettings("KEY_SPOT_COLORING", color);

		displayer.refresh();
		displayer.render();
	}
}
