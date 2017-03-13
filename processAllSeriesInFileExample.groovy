
// @ImageJ ij

import java.util.HashMap;
import java.util.concurrent.Future;
import uk.ac.cruk.SatelliteAnalysis;
import ij.ImagePlus;
import ij.IJ;
import ij.ImageStack;
import ij.process.ImageProcessor;
import ij.measure.Calibration;
import ij.plugin.HyperStackConverter;
import loci.formats.ChannelSeparator;
import loci.formats.meta.MetadataRetrieve;
import loci.formats.services.OMEXMLServiceImpl;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.util.LociPrefs;
import ome.units.UNITS;

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


String id = "C:\\Users\\pike01\\Documents\\Side_Projects\\Data\\Valentina_Quarantotti\\PCM1\\RPE1_siRNA exp12-09-2016_PCM1-488_gTub-555_28-09-2016\\RPE1_siRNA exp12-09-2016_PCM1-488_gTub-555_28-09-2016.lif";

ImageProcessorReader r = new ImageProcessorReader(new ChannelSeparator(LociPrefs.makeImageReader()));
OMEXMLServiceImpl OMEXMLService = new OMEXMLServiceImpl();
r.setMetadataStore(OMEXMLService.createOMEXMLMetadata());
r.setId(id);
MetadataRetrieve meta = (MetadataRetrieve) r.getMetadataStore();

for (int i = 0; i < 30; i = i + 10) {
	println(i);
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
	fc = ij.command().run(SatelliteAnalysis.class, true, pluginSettings);
	while (!fc.isDone()){
	}
}
r.close();



