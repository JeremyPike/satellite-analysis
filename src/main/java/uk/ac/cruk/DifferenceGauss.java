

package uk.ac.cruk;

import org.scijava.ItemIO;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;
import org.scijava.ui.UIService;
import net.imagej.Dataset;
import net.imagej.DatasetService;
import net.imagej.ImageJ;
import net.imagej.ops.OpService;
import net.imglib2.IterableInterval;
import net.imglib2.img.Img;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.real.FloatType;

/**
 * This is a minimal ImageJ command implementing a difference of Gaussian filter.
 */

@Plugin(type = Command.class, menuPath = "Plugins>DoG Filtering")
public class DifferenceGauss<T extends RealType<T>> implements Command {
	
	
	@Parameter
	private double smallSigma;
	
	@Parameter
	private double bigSigma;
	
	@Parameter
	private Dataset currentData;

	@Parameter
	private UIService uiService;

	@Parameter
	private OpService opService;

	@Parameter
	private DatasetService datasetService;

	@Parameter(type = ItemIO.OUTPUT)
	private IterableInterval<FloatType> dog;

	@Override
	public void run() {
		final Img<T> image = (Img<T>) currentData.getImgPlus();

		// Convert image to FloatType for better numeric precision
		Img<FloatType> converted = opService.convert().float32(image);

		// Create two images and peform Gaussian filtering
		Img<FloatType> smallGauss = opService.create().img(converted);
		Img<FloatType> bigGauss = opService.create().img(converted);
		opService.filter().gauss(smallGauss, converted, smallSigma);
		opService.filter().gauss(bigGauss, converted, bigSigma);

		// uiService.show(smallGauss);
		// uiService.show(bigGauss);
		
		// calculate difference of Gaussian output
		dog = opService.create().img(converted);
		opService.run("math.subtract", dog, smallGauss, bigGauss);


	}

	/**
	 * This main function serves for development purposes.
	 *
	 * @param args
	 *            whatever, it's ignored
	 * @throws Exception
	 */
	public static void main(final String... args) throws Exception {
		
		// create the ImageJ application context
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// load example dataset
		final Dataset dataset = ij.scifio().datasetIO().open("http://imagej.net/images/FluorescentCells.jpg");
		
		// show the image
		ij.ui().show(dataset);

		// invoke the plugin
		ij.command().run(DifferenceGauss.class, true);

	}

}
