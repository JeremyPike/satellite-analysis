
// @ImageJ ij
// @File(label="Select a file") file

import java.util.HashMap;
import java.util.concurrent.Future;
import uk.ac.cruk.SatelliteAnalysis;
import ij.ImagePlus;
import ij.IJ;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.HyperStackConverter;
import loci.formats.ChannelSeparator;
import loci.formats.meta.MetadataRetrieve;
import loci.formats.services.OMEXMLServiceImpl;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.util.LociPrefs;
import ome.units.UNITS;
import org.scijava.command.CommandModule;

// Create emtpy map for plugin settings
Map<String, Object> pluginSettings = new HashMap<String, Object>();

// Set plugin settings
// Radius for Gaussian blur (DAPI channel) (microns)
pluginSettings.put("dapiSigmaMicrons", 1.0d);
// Radius for nuclei spot detection (microns)
pluginSettings.put("dapiSpotRadius", 7.0d);
// Quality threshold for nuclei spots
pluginSettings.put("dapiSpotQualityThresh", 0.075d);
// Radius for satellite spot detection (microns
pluginSettings.put("satSpotRadius", 0.2d);
// Mean intenisty threshold for satelllite spots
pluginSettings.put("satSpotMeanIntThresh", 26.0d);
// Radius for Gaussian blur (Centrosome channel) (microns)
pluginSettings.put("centSigmaMicrons", 0.2d);
// Radius for rolling ball BS (Centrosome channel) (microns)
pluginSettings.put("centRollingBallRad", 2.0d);
// Radius for centrosome spot detection (microns)
pluginSettings.put("centSpotRadius", 0.75d);
// Max intenisty threshold for centrosome spots
pluginSettings.put("centSpotMaxIntThresh", 18.0d);
// Division factor for measurement within DAPI spots
pluginSettings.put("dapiCentreFrac", 4.0d);
// Multiplication edge factor for minimum distance from image edge for Dapi spots
pluginSettings.put("edgeFac", 2.0d);
// Multiplication factor for median to calculate maximum nuceli medan intensity
pluginSettings.put("nuceliMedianFac", 1.5d);
// Filter satellites closer to centrosomes than centrosome spot radius
pluginSettings.put("filterSatellites", true);
// Display spot detection overlays on maximal projection
pluginSettings.put("displayMPOverlay", false);

// Reader to load specified file
ImageProcessorReader r = new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
OMEXMLServiceImpl OMEXMLService = new OMEXMLServiceImpl();
r.setMetadataStore(OMEXMLService.createOMEXMLMetadata());
r.setId(file.getPath());

// get metadata
MetadataRetrieve meta = (MetadataRetrieve) r.getMetadataStore();

// get current results table and reset
ResultsTable rt = ResultsTable.getResultsTable();
rt.reset();

// get number of series in specified file
int numSeries = meta.getImageCount();

// iterate through all series
for (int i = 0; i < numSeries; i++) {
	println("Processing series " + (i + 1) + " of " + numSeries);
	// set reader to current series
	r.setSeries(i);
	// load data for current series
	ImageStack stack = new ImageStack(r.getSizeX(), r.getSizeY());
	for (int n = 0; n < r.getImageCount(); n++) {
		ImageProcessor ip = r.openProcessors(n)[0];
		stack.addSlice("" + (n + 1), ip);
	}
	ImagePlus imp = new ImagePlus("", stack);
	// check dimension are correct
	imp = HyperStackConverter.toHyperStack(imp, r.getSizeC(), r.getSizeZ(), r.getSizeT());
	// set scale using metadata
	Calibration cali = new Calibration();
	cali.pixelWidth = meta.getPixelsPhysicalSizeX(i).value(UNITS.MICROMETER).doubleValue();
	cali.pixelHeight = meta.getPixelsPhysicalSizeY(i).value(UNITS.MICROMETER).doubleValue();
	cali.pixelDepth = meta.getPixelsPhysicalSizeZ(i).value(UNITS.MICROMETER).doubleValue();
	cali.setUnit("micron");
	imp.setGlobalCalibration(cali);
	// set window title to series name
	imp.setTitle(meta.getImageName(i));
	// add the series data to the plugin settings
	pluginSettings.put("currentData", imp);
	// run the satellite analysis command
	Future<CommandModule> fc = ij.command().run(SatelliteAnalysis.class, true, pluginSettings);
	// wait for command to finish before continuing 
	while (!fc.isDone()){
	}
}
// close the reader
r.close();



