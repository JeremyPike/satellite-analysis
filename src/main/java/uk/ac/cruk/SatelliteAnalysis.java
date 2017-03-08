
package uk.ac.cruk;

import java.util.Map;

import org.scijava.command.Command;
import org.scijava.display.Display;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import fiji.plugin.trackmate.Model;
import fiji.plugin.trackmate.SelectionModel;
import fiji.plugin.trackmate.Settings;
import fiji.plugin.trackmate.SpotCollection;
import fiji.plugin.trackmate.TrackMate;
import fiji.plugin.trackmate.detection.DownsampleLogDetectorFactory;
import fiji.plugin.trackmate.features.FeatureFilter;
import fiji.plugin.trackmate.visualization.hyperstack.HyperStackDisplayer;
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
	
	@Parameter(label = "Quality threshold for nuclei spot quality: ", persist = false, min = "0")
	private double dapiSpotQualityThresh = 0.075;
	
	
	
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
		//Display<?> dapiImageMPBlurDisplay = ij.display().createDisplay(dapiImageMPBlur);

		ImagePlus dapiImageMPBlurImp = legacy.getImageMap().registerDataset(datasetService.create(dapiImageMPBlur));
		
		// Use Trackmate to count the number of nuclei
		Settings settings = new Settings();
		settings.detectorFactory = new DownsampleLogDetectorFactory();
		settings.addSpotFilter(new FeatureFilter("QUALITY", dapiSpotQualityThresh, true));
		Map<String, Object> map = settings.detectorFactory.getDefaultSettings();
		map.put("RADIUS", dapiSpotRadius);
		map.put("DO_MEDIAN_FILTERING", false);
		map.put("DOWNSAMPLE_FACTOR", 4);
		map.put("DO_SUBPIXEL_LOCALIZATION", false);
		map.put("TARGET_CHANNEL", 1.0d);
		map.put("THRESHOLD", 0.0d);
		settings.detectorSettings = map;
		settings.setFrom(dapiImageMPBlurImp);
        
		SpotCollection dapiSpots = countSpotsTrackmate(settings, false);
	
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
