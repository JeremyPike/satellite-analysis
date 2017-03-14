
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


Map<String, Object> pluginSettings = new HashMap<String, Object>();
 
pluginSettings.put("dapiSigmaMicrons", 1.0d);
pluginSettings.put("dapiSpotRadius", 7.0d);
pluginSettings.put("dapiSpotQualityThresh", 0.075d);
pluginSettings.put("satSpotRadius", 0.2d);
pluginSettings.put("satSpotMeanIntThresh", 26.0d);
pluginSettings.put("centSigmaMicrons", 0.2d);
pluginSettings.put("centRollingBallRad", 2.0d);
pluginSettings.put("centSpotRadius", 0.75d);
pluginSettings.put("centSpotMaxIntThresh", 18.0d);
pluginSettings.put("dapiCentreFrac", 4.0d);
pluginSettings.put("edgeFac", 2.0d);
pluginSettings.put("nuceliMedianFac", 1.5d);
pluginSettings.put("filterSatellites", true);
pluginSettings.put("displayMPOverlay", false);


ImageProcessorReader r = new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
OMEXMLServiceImpl OMEXMLService = new OMEXMLServiceImpl();
r.setMetadataStore(OMEXMLService.createOMEXMLMetadata());
r.setId(file.getPath());
MetadataRetrieve meta = (MetadataRetrieve) r.getMetadataStore();

ResultsTable rt = ResultsTable.getResultsTable();
rt.reset();

int numSeries = meta.getImageCount();

for (int i = 0; i < numSeries; i = i + 1) {
	println("Processing series " + (i + 1) + " of " + numSeries);
	r.setSeries(i);
	ImageStack stack = new ImageStack(r.getSizeX(), r.getSizeY());
	for (int n = 0; n < r.getImageCount(); n++) {
		ImageProcessor ip = r.openProcessors(n)[0];
		stack.addSlice("" + (n + 1), ip);
	}
	ImagePlus imp = new ImagePlus("", stack);
	imp = HyperStackConverter.toHyperStack(imp, r.getSizeC(), r.getSizeZ(), r.getSizeT());
	
	Calibration cali = new Calibration();
	cali.pixelWidth = meta.getPixelsPhysicalSizeX(i).value(UNITS.MICROMETER).doubleValue();
	cali.pixelHeight = meta.getPixelsPhysicalSizeY(i).value(UNITS.MICROMETER).doubleValue();
	cali.pixelDepth = meta.getPixelsPhysicalSizeZ(i).value(UNITS.MICROMETER).doubleValue();
	cali.setUnit("micron");
	imp.setGlobalCalibration(cali);
	
	imp.setTitle(meta.getImageName(i));
	pluginSettings.put("currentData", imp);
	Future<CommandModule> fc = ij.command().run(SatelliteAnalysis.class, true, pluginSettings);
	while (!fc.isDone()){
	}
}
r.close();



