
package uk.ac.cruk;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.HyperStackConverter;
import ij.process.ImageProcessor;
import loci.formats.ChannelSeparator;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.util.LociPrefs;
import net.imagej.Data;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.axis.Axes;
import net.imagej.legacy.LegacyService;
import net.imagej.ops.Op;
import net.imagej.ops.OpService;
import net.imagej.ops.Ops;
import net.imagej.ops.special.computer.UnaryComputerOp;
import net.imglib2.Dimensions;
import net.imglib2.FinalInterval;
import net.imglib2.IterableInterval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.DoubleType;
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
	private UIService uiService;

	@Parameter
	private OpService opService;

	@Parameter
	private DatasetService datasetService;

	@Parameter
	private double sigma;

	@Override
	public void run() {

		final Img<T> image = (Img<T>) currentData.getImgPlus();

		long width = image.dimension(0);
		long height = image.dimension(1);
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
		opService.filter().gauss(dapiImageMPBlur, dapiImageMPBlur, sigma);
		uiService.show(dapiImageMPBlur);

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
		r.setId(id);
		r.setSeries(0);
		ImageStack stack = new ImageStack(r.getSizeX(), r.getSizeY());
		for (int n = 0; n < r.getImageCount(); n++) {
			ImageProcessor ip = r.openProcessors(n)[0];
			stack.addSlice("" + (n + 1), ip);
		}
		ImagePlus imp = new ImagePlus("", stack);
		imp = HyperStackConverter.toHyperStack(imp, r.getSizeC(), r.getSizeZ(), r.getSizeT());
		imp.show();
		r.close();
		// invoke the plugin
		ij.command().run(SatelliteAnalysis.class, true);

	}

}
