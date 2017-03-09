
package uk.ac.cruk;

import java.util.Iterator;
import java.util.Map;

import org.scijava.command.Command;
import org.scijava.display.Display;
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
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer;
import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.measure.Calibration;
import ij.plugin.HyperStackConverter;
import ij.process.ImageProcessor;
import loci.formats.ChannelSeparator;
import loci.formats.meta.MetadataRetrieve;
import loci.formats.services.OMEXMLServiceImpl;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.util.LociPrefs;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.display.ImageDisplay;
import net.imagej.legacy.LegacyService;
import net.imagej.ops.Op;
import net.imagej.ops.OpService;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;
import ome.units.UNITS;

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
	
	@Override
	public void run() {
		// Create legacy service
		LegacyService legacy = ij.get(LegacyService.class);
		
		// Get IJ2 Img object
		final Img<T> image = (Img<T>) currentData.getImgPlus();

		// Pull out data size
		long width = image.dimension(0);
		long height = image.dimension(1);
		long numSlices = image.dimension(3);


		// Get Dapi data
		FinalInterval slice = FinalInterval.createMinSize(0, 0, 0, 0, width, height, 1, numSlices);
		RandomAccessibleInterval<T> dapiImage = opService.transform().crop(image, slice);
		
		// Dimensions for max proj image
		long[] projectedDimensions = new long[] {width, height};
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
		IJ.run(dapiImageMPBlurImp, "Properties...", "channels=1 slices=1 frames=1 unit=micron pixel_width=" +  currentData.averageScale(0) + " pixel_height=" +  currentData.averageScale(1) + " voxel_depth=1");

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
		IJ.run(satImageImp, "Properties...", "channels=1 slices=" + image.dimension(3) + " frames=" + image.dimension(4) + " unit=micron pixel_width=" +  currentData.averageScale(0) + " pixel_height=" +  currentData.averageScale(1) + " voxel_depth=" +  currentData.averageScale(3));
	
		
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
		
	//	Img<FloatType> centImageBlur = opService.convert().float32((IterableInterval<C>) centImage);

		// Gaussian filtering
		opService.filter().gauss(centImageBlur, centImageBlur, centSigmaMicrons / currentData.averageScale(0), centSigmaMicrons / currentData.averageScale(1), centSigmaMicrons / currentData.averageScale(3));
		
		// Convert to IJ1 ImagePlus for Trackmate
		ImagePlus centImageBlurImp = legacy.getImageMap().registerDataset(datasetService.create(centImageBlur));
		IJ.run(centImageBlurImp, "Properties...", "channels=1 slices=" + image.dimension(3) + " frames=" + image.dimension(4) + " unit=micron pixel_width=" +  currentData.averageScale(0) + " pixel_height=" +  currentData.averageScale(1) + " voxel_depth=" +  currentData.averageScale(3));
		// Rolling ball background subtraction
		IJ.run(centImageBlurImp, "Subtract Background...", "rolling=" + 20 + " stack");
		
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
        
		SpotCollection centSpots = countSpotsTrackmate(centSettings, false);

		int numNuclei = dapiSpots.getNSpots(true);
		Iterator<Spot> dapiIterator = dapiSpots.iterator(true);
		Nucleus[] nuclei = new Nucleus[numNuclei];
		for (int i = 0; i < numNuclei; i++) {
			Spot dapiSpot = dapiIterator.next();
			nuclei[i] = new Nucleus(dapiSpot);
		}
		
		Iterator<Spot> centIterator = centSpots.iterator(true);
		Spot cent;	
		while (centIterator.hasNext()) {
			cent = centIterator.next();
			int closestID = dapiSpots.getClosestSpot(cent, 0, true).ID();
			for (int i = 0; i < numNuclei; i++) {
				if (closestID == nuclei[i].getNucleusSpot().ID()) {
					nuclei[i].addCentrosome(cent);
				}
			}
		}
		Iterator<Spot> satIterator = satSpots.iterator(true);
		Spot sat;	
		while (satIterator.hasNext()) {
			sat = satIterator.next();
			int closestID = dapiSpots.getClosestSpot(sat, 0, true).ID();
			for (int i = 0; i < numNuclei; i++) {
				if (closestID == nuclei[i].getNucleusSpot().ID()) {
					nuclei[i].addSatellite(sat);
				}
			}
		}
		
		nuclei[2].display(centImageBlurImp);
		
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
		String id = "C:\\Users\\pike01\\Documents\\testSatSeries.lif";
		ImageProcessorReader r = new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
		OMEXMLServiceImpl OMEXMLService = new OMEXMLServiceImpl();
		r.setMetadataStore(OMEXMLService.createOMEXMLMetadata());
		r.setId(id);
		r.setSeries(0);
		ImageStack stack = new ImageStack(r.getSizeX(), r.getSizeY());
		for (int n = 0; n < r.getImageCount(); n++) {
			ImageProcessor ip = r.openProcessors(n)[0];
			stack.addSlice("" + (n + 1), ip);
		}
		ImagePlus imp = new ImagePlus("", stack);
		imp = HyperStackConverter.toHyperStack(imp, r.getSizeC(), r.getSizeZ(), r.getSizeT());
		MetadataRetrieve meta = (MetadataRetrieve) r.getMetadataStore();
		Calibration cali = new Calibration();
		cali.pixelWidth = meta.getPixelsPhysicalSizeX(0).value(UNITS.MICROMETER).doubleValue();
		cali.pixelHeight = meta.getPixelsPhysicalSizeY(0).value(UNITS.MICROMETER).doubleValue();
		cali.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value(UNITS.MICROMETER).doubleValue();
		cali.setUnit("micron");
		imp.setGlobalCalibration(cali);
		imp.show();
		r.close();
		// invoke the plugin
		ij.command().run(SatelliteAnalysis.class, true);

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
			HyperStackDisplayer displayer =  new HyperStackDisplayer(model, selectionModel, settings.imp);
			displayer.render();
			displayer.refresh();
	    }
		//Return spot collection
	    return model.getSpots();
	}
	
	
}
