
package uk.ac.cruk;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.Future;

import org.python.antlr.PythonParser.printlist2_return;
import org.scijava.command.Command;
import org.scijava.command.CommandModule;
import org.scijava.module.ModuleItem;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.DialogPrompt.MessageType;
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
 * The TrackMate Fiji plugin is used for spot detection: Tinevez, Jean-Yves, et
 * al.
 * "TrackMate: An open and extensible platform for single-particle tracking."
 * Methods (2016).
 * 
 * @author Jeremy Pike
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

		// create legacy service for conversion between IJ and IJ2 classes
		LegacyService legacy = ij.get(LegacyService.class);

		// get the name of the current data
		String seriesName = currentData.getName();

		// get IJ2 Img object
		final Img<T> image = (Img<T>) currentData.getImgPlus();

		// pull out data dimensions
		long width = image.dimension(0);
		long height = image.dimension(1);
		long numChannels = image.dimension(2);
		long numSlices = image.dimension(3);

		// check dataset has three channels
		if (numChannels != 3) {
			uiService.showDialog("Please privide data with three channels (DAPI, satellite, centrosome)",
					MessageType.ERROR_MESSAGE);
			return;
		}

		// crop Dapi channel
		FinalInterval slice = FinalInterval.createMinSize(0, 0, 0, 0, width, height, 1, numSlices);
		RandomAccessibleInterval<T> dapiImage = opService.transform().crop(image, slice);

		// dimensions for max proj image
		long[] projectedDimensions = new long[] { width, height };
		// generate memory for max proj with type same as input data
		Img<T> dapiImageMP = (Img<T>) opService.run("create.img", projectedDimensions, image.firstElement());
		// maximal projection along z axis
		Op maxProj = opService.op("stats.max", dapiImage);
		opService.run("transform.project", dapiImageMP, dapiImage, opService.op("stats.max", dapiImage), 2);

		// convert image to FloatType for better numeric precision
		Img<FloatType> dapiImageMPBlur = opService.convert().float32(dapiImageMP);
		// Gaussian filtering
		opService.filter().gauss(dapiImageMPBlur, dapiImageMPBlur, dapiSigmaMicrons / currentData.averageScale(0));

		// convert to IJ1 ImagePlus objects for Trackmate
		// for spot detection
		ImagePlus dapiImageMPBlurImp = legacy.getImageMap().registerDataset(datasetService.create(dapiImageMPBlur));
		// raw data for mean DAPI measurements
		ImagePlus dapiImageMPImp = legacy.getImageMap().registerDataset(datasetService.create(dapiImageMP));
		// check calibration and dimensions are in agreement with original
		// dataset
		IJ.run(dapiImageMPBlurImp, "Properties...", "channels=1 slices=1 frames=1 unit=micron pixel_width="
				+ currentData.averageScale(0) + " pixel_height=" + currentData.averageScale(1) + " voxel_depth=1");
		IJ.run(dapiImageMPImp, "Properties...", "channels=1 slices=1 frames=1 unit=micron pixel_width="
				+ currentData.averageScale(0) + " pixel_height=" + currentData.averageScale(1) + " voxel_depth=1");

		// set up trackmate for dapi channel
		Settings dapiDettings = new Settings();
		// downsample Log detector
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

		// use trackmate to detect DAPI spots
		SpotCollection dapiSpots = countSpotsTrackmate(dapiDettings, false);

		// crop satellite data
		slice = FinalInterval.createMinSize(0, 0, 1, 0, width, height, 1, numSlices);
		RandomAccessibleInterval<T> satImage = opService.transform().crop(image, slice);
		// Convert image to FloatType for better numeric precision
		Img<FloatType> satImageFloat = (Img<FloatType>) opService.run("convert.float32", satImage);

		// convert to IJ1 ImagePlus objects for Trackmate
		ImagePlus satImageImp = legacy.getImageMap().registerDataset(datasetService.create(satImageFloat));
		// check calibration and dimensions are in agreement with original
		// dataset
		IJ.run(satImageImp, "Properties...",
				"channels=1 slices=" + image.dimension(3) + " frames=" + image.dimension(4)
						+ " unit=micron pixel_width=" + currentData.averageScale(0) + " pixel_height="
						+ currentData.averageScale(1) + " voxel_depth=" + currentData.averageScale(3));

		// Use Trackmate to count the number of satellites
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

		// crop centrosome data
		slice = FinalInterval.createMinSize(0, 0, 2, 0, width, height, 1, numSlices);
		RandomAccessibleInterval<T> centImage = opService.transform().crop(image, slice);
		// convert image to FloatType for better numeric precision
		Img<FloatType> centImageBlur = (Img<FloatType>) opService.run("convert.float32", centImage);

		// Gaussian filtering
		opService.filter().gauss(centImageBlur, centImageBlur, centSigmaMicrons / currentData.averageScale(0),
				centSigmaMicrons / currentData.averageScale(1), centSigmaMicrons / currentData.averageScale(3));

		// convert to IJ1 ImagePlus for Trackmate
		// for spot detection
		ImagePlus centImageBlurImp = legacy.getImageMap().registerDataset(datasetService.create(centImageBlur));
		// raw data for mean centrosome intensity measurements
		ImagePlus centImageImp = legacy.getImageMap().registerDataset(datasetService.create(centImage));
		// check calibration and dimensions are in agreement with original
		// dataset
		IJ.run(centImageBlurImp, "Properties...",
				"channels=1 slices=" + image.dimension(3) + " frames=" + image.dimension(4)
						+ " unit=micron pixel_width=" + currentData.averageScale(0) + " pixel_height="
						+ currentData.averageScale(1) + " voxel_depth=" + currentData.averageScale(3));
		IJ.run(centImageImp, "Properties...",
				"channels=1 slices=" + image.dimension(3) + " frames=" + image.dimension(4)
						+ " unit=micron pixel_width=" + currentData.averageScale(0) + " pixel_height="
						+ currentData.averageScale(1) + " voxel_depth=" + currentData.averageScale(3));
		// rolling ball background subtraction
		IJ.run(centImageBlurImp, "Subtract Background...",
				"rolling=" + centRollingBallRad / currentData.averageScale(0) + " stack");

		// Use Trackmate to count the number of centrosomes
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

		// for each detected dapi spots generate a Nucelus object and pull out
		// the mean dapi itensity
		int numUnfilteredNuclei = dapiSpots.getNSpots(true);
		Iterator<Spot> dapiIterator = dapiSpots.iterator(true);
		Nucleus[] nucleiUnfiltered = new Nucleus[numUnfilteredNuclei];
		double[] nuceliMean = new double[numUnfilteredNuclei];
		for (int i = 0; i < numUnfilteredNuclei; i++) {
			Spot dapiSpot = dapiIterator.next();
			nucleiUnfiltered[i] = new Nucleus(dapiSpot, dapiImageMPImp, dapiCentreFrac, edgeFac);
			nuceliMean[i] = nucleiUnfiltered[i].getMeanDapiIntensity();
		}

		// calculate the median mean dapit intensity. This is used later to
		// identify mitotic cells
		double nucleiMeanMedian;
		Arrays.sort(nuceliMean);
		if (nuceliMean.length % 2 == 0) {
			nucleiMeanMedian = (nuceliMean[nuceliMean.length / 2] + nuceliMean[nuceliMean.length / 2 - 1]) / 2;
		} else {
			nucleiMeanMedian = nuceliMean[nuceliMean.length / 2];
		}

		// for each dapi detection find the distance to every other detection.
		// Only include in the filtered nuclei ArrayList if it has higher mean
		// intensity than all detections closer than a threshold distance

		// array list for filtered Nucelus objects
		ArrayList<Nucleus> nucleiAL = new ArrayList<Nucleus>();

		// loop through all unfiltered nuclei
		for (int i = 0; i < numUnfilteredNuclei; i++) {

			boolean nucCheck = true;
			// find distance to all other nuceli
			for (int j = 0; j < numUnfilteredNuclei; j++) {
				double nucSepDist = Nucleus.spotSeperation(nucleiUnfiltered[i].getNucleusSpot(),
						nucleiUnfiltered[j].getNucleusSpot(), true);
				// if closer than user defined threshold
				if (i != j && nucSepDist <= dapiSpotRadius
						&& nucleiUnfiltered[i].getMeanDapiIntensity() < nucleiUnfiltered[j].getMeanDapiIntensity()) {
					// set check to false
					nucCheck = false;
					// change type for colour displays
					nucleiUnfiltered[i].setNucleusSpotType(1.0d);
				}

			}
			// if check is true include in filtered ArrayList
			if (nucCheck) {
				nucleiAL.add(nucleiUnfiltered[i]);
			}
		}

		// convert ArrayList to array and find numbered of filtered nuceli
		Nucleus[] nuclei = nucleiAL.toArray(new Nucleus[nucleiAL.size()]);
		int numNuclei = nuclei.length;

		// assign each detected centrosome to a filtered Nucelus based on
		// closest distance
		Iterator<Spot> centIterator = centSpots.iterator(true);
		Spot cent;
		double minDist;
		int minInd;
		// iterate through all centrosome Spots
		while (centIterator.hasNext()) {
			cent = centIterator.next();
			// find minimum distance (2D) and corresponding array index
			minDist = Nucleus.spotSeperation(cent, nuclei[0].getNucleusSpot(), true);
			minInd = 0;
			for (int i = 1; i < numNuclei; i++) {
				double dist = Nucleus.spotSeperation(cent, nuclei[i].getNucleusSpot(), true);
				if (dist < minDist) {
					minDist = dist;
					minInd = i;
				}
			}
			// add centrosome to correct Nucleus
			nuclei[minInd].addCentrosome(cent);
		}
	
		// assign each detected satellite to a filtered Nucelus based on closest
		// distance
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

		// get current ResultsTable
		ResultsTable rt = ResultsTable.getResultsTable();

		// loop through all filtered nuceli
		for (int i = 0; i < numNuclei; i++) {

			// if user has specified satellites close to centrosomes should be
			// removed from the analysis
			if (filterSatellites) {
				// compute Nucleus properties and statistics
				nuclei[i].computeProperties(centSpotRadius, nucleiMeanMedian * nuceliMedianFac);
			} else {
				nuclei[i].computeProperties(0.0d, nucleiMeanMedian * nuceliMedianFac);
			}
			// if Nucleus isnt close to image edges and isnt labelled as mitotic
			// add its properties to the ResultTable
			if (!nuclei[i].getIsEdge() && !nuclei[i].getIsMitotic()) {
				rt.incrementCounter();
				// set spot name for nucleus based on current results table row
				nuclei[i].setNucleusSpotName(Integer.toString(rt.getCounter()));
				// include series name
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
		// show and update ResultsTable
		rt.show("Results");

		// if user has specified a detection visualisation should be generated
		if (displayMPOverlay) {

			// dimensions for max proj image of original dataset
			projectedDimensions = new long[] { width, height, numChannels };
			// generate memory for max proj with type same as input data
			Img<T> allImageMP = (Img<T>) opService.run("create.img", projectedDimensions, image.firstElement());
			// maximal projection along z axis
			maxProj = opService.op("stats.max", currentData);
			opService.run("transform.project", allImageMP, currentData, opService.op("stats.max", currentData), 3);
			// convert to IJ1 ImagePlus
			ImagePlus allImageMPImp = legacy.getImageMap().registerDataset(datasetService.create(allImageMP));
			// check calibration and dimensions are in agreement with original
			// dataset
			IJ.run(allImageMPImp, "Properties...",
					"channels=" + numChannels + " slices=1 frames=1 unit=micron pixel_width="
							+ currentData.averageScale(0) + " pixel_height=" + currentData.averageScale(1)
							+ " voxel_depth=1");
			// create CompositeImage so LUTs can be set and shown as composite
			// image
			CompositeImage impComp = new CompositeImage(allImageMPImp);
			impComp.setChannelLut(LUT.createLutFromColor(Color.BLUE), 1);
			impComp.setChannelLut(LUT.createLutFromColor(Color.GREEN), 2);
			impComp.setChannelLut(LUT.createLutFromColor(Color.RED), 3);
			allImageMPImp = (ImagePlus) impComp.clone();
			allImageMPImp.setDisplayMode(IJ.COMPOSITE);

			// create detection visualisation
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

		// create imageJ and show UI
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// open test dataset and show
		String id = "C:\\Users\\pike01\\Documents\\testSatSeries.tif";
		Dataset dataset = (Dataset) ij.io().open(id);
		ij.ui().show(dataset);

		// get default setting input parameters for command and save in Map
		Map<String, Object> pluginSettings = new HashMap<String, Object>();
		Iterator<ModuleItem<?>> inputs = ij.command().getCommand(SatelliteAnalysis.class).inputs().iterator();
		while (inputs.hasNext()) {
			ModuleItem<?> input = inputs.next();
			if (input.getDefaultValue() != null)
				pluginSettings.put(input.getName(), input.getDefaultValue());
		}
		// run command with default settings
		Future<CommandModule> fc = ij.command().run(SatelliteAnalysis.class, true, pluginSettings);

	}

	/**
	 * Uses trackmate to detect spots
	 *
	 * @param settings
	 *            Settings for trackmate spot detection
	 * 
	 * @param display
	 *            if true displays the results
	 * 
	 * @return The collection of spots from trackmate
	 */
	public static SpotCollection countSpotsTrackmate(Settings settings, boolean display) {

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
		// display spot detections
		if (display) {
			SelectionModel selectionModel = new SelectionModel(model);
			HyperStackDisplayer displayer = new HyperStackDisplayer(model, selectionModel, settings.imp);
			displayer.render();
			displayer.refresh();
		}
		// Return spot collection
		return model.getSpots();
	}

	/**
	 * Uses trackmate to detect spots
	 *
	 * @param settings
	 *            Settings for trackmate spot detection
	 * 
	 * @param display
	 *            if true displays the results
	 * 
	 * @param imp
	 *            uses this image for computation of features after spot
	 *            detection
	 * 
	 * @return The collection of spots from trackmate
	 */
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

	/**
	 * displays all nuclear, centrosome and satellite detections contained in
	 * the given array of Nucleus objects
	 *
	 * @param imp
	 *            ImagePlus used to display detections. This should only have
	 *            one slice
	 * 
	 * @param nuclei
	 *            array of nucleus objects containing the relevant detections
	 * 
	 * @param filtered
	 *            if true will only show filtered satellites
	 * 
	 */
	private void displayAll(ImagePlus imp, Nucleus[] nuclei, boolean filtered) {

		// create spot collection to hold all detections
		SpotCollection allSpots = new SpotCollection();
		Spot spot;
		SpotCollection satellites;
		SpotCollection centrosomes;
		Iterator<Spot> satIterator;
		Iterator<Spot> centIterator;

		// loop through all Nucleus objects
		for (int i = 0; i < nuclei.length; i++) {

			// got nuclear spot
			spot = nuclei[i].getNucleusSpot();
			// set z position to zero as visualisation is 2D
			spot.putFeature("POSITION_Z", 0.0d);
			// add to joint collection
			allSpots.add(spot, 0);
			// either get the filtered satellites or the unfiltered
			if (filtered) {
				satellites = nuclei[i].getSatellitesFilt();
			} else {
				satellites = nuclei[i].getSatellites();
			}
			centrosomes = nuclei[i].getCentrosomes();
			// add all centrosomes and satellites in given Nucleus to joint
			// collection
			satIterator = satellites.iterator(true);
			while (satIterator.hasNext()) {
				spot = satIterator.next();
				// set z position to zero as visualisation is 2D
				spot.putFeature("POSITION_Z", 0.0d);
				allSpots.add(spot, 0);
			}
			centIterator = centrosomes.iterator(true);
			while (centIterator.hasNext()) {
				spot = centIterator.next();
				// set z position to zero as visualisation is 2D
				spot.putFeature("POSITION_Z", 0.0d);
				allSpots.add(spot, 0);
			}

		}
		
	
		// Add joint spot collection to an empty Trackmate model
		Model modelAll = new Model();
		modelAll.setSpots(allSpots, false);

		// Create spot colour scheme based on "TYPE" property
		SpotColorGenerator color = new SpotColorGenerator(modelAll);
		color.setFeature("TYPE");
		// TYPE can range from 0 to 5 for all spot types
		color.setMinMax(0, 5);
		// Display spots on provided ImagePlus using a Hyperstack displayer
		HyperStackDisplayer displayer = new HyperStackDisplayer(modelAll, new SelectionModel(modelAll), imp);
		displayer.setDisplaySettings(HyperStackDisplayer.KEY_SPOT_COLORING, color);
		// display spot names which correspond to the results table row
		displayer.setDisplaySettings(HyperStackDisplayer.KEY_DISPLAY_SPOT_NAMES, true);
		displayer.refresh();
		displayer.render();
	}
}
