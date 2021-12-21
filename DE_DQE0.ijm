print("\\Clear");
print("DE_DQE0");
print("Version : 1.01");
print("Date : 2021.12.21");
print("Author : Benjamin Bammes, Direct Electron LP (bbammes@directelectron.com)");
print("License : GNU General Public License v2.0");
print("");
print("Uses the noise binning method as described in McMullan, et al. 2009. Ultramicroscopy 109:1126-43.");
print("");

openImages = getList("window.titles");
if (openImages.length < 1) {
	print("Error");
	exit("No image stack found.");
}
Stack.getDimensions(width, height, channels, slices, frames);
if (slices < 2) {
	print("Error");
	exit("This macro requires an image stack with at least two slices.");
}
if (getValue("selection.size") != 4) {
	print("Error");
	exit("You must select a single square ROI before running this macro.");
}

run("Duplicate...", "title=DQE0_CroppedStack duplicate");
Stack.getDimensions(width, height, channels, slices, frames);
if (width != height) {
	close("DQE0_CroppedStack");
	print("Error");
	exit("You must select a square ROI before running this macro.");
}

aduPerElectron = getNumber("ADU per Primary Electron", 0);
if (aduPerElectron < 1.0) {
	close("DQE0_CroppedStack");
	print("Error");
	exit("You must enter a valid value (>=1) for ADU per Primary Electron.");
}

maxBinning = Math.floor(width / 14);

resultArrayX = newArray(0);
resultArrayY = newArray(0);
maxResultValue = 0.0;

Stack.getStatistics(voxelCount, mean, min, max, std);
minValue = mean - std * 4.0;
maxValue = mean + std * 4.0;
run("Min...", "value=" + minValue + " stack");
run("Max...", "value=" + maxValue + " stack");
Stack.getStatistics(voxelCount, meanADU, min, max, std);
electronsPerPixel = meanADU / aduPerElectron;

for (i = 1; i <= slices; i++) {
	Stack.setSlice(i);
	run("Set Label...", "label=DQE0_Image_" + i);
}
run("Stack to Images");

Plot.create("DQE0 Noise Binning Results", "Binning", "Variance / Binning^2");
Plot.setFrameSize(800, 600);
Plot.setLimits(0, maxBinning, 0, 100 );
Plot.setColor("#1965B0");
Plot.add("circles", resultArrayX, resultArrayY);
Plot.show;

comparisons = 0;

for (i = 1; i <= slices; i++) {
	for (j = i + 1; j <= slices; j++) {
		print("Processing frame " + i + " and " + j + "...");
		imageCalculator("Subtract create 32-bit", "DQE0_Image_" + i, "DQE0_Image_" + j);
		selectWindow("Result of DQE0_Image_" + i);
		rename("DQE0_ActiveImage");
		for (b = 1; b < maxBinning; b++) {
			selectWindow("DQE0_ActiveImage");
			run("Duplicate...", "title=DQE0_ActiveImage_Binned");
			run("Bin...", "x=" + b + " y=" + b + " bin=Sum");
			getDimensions(binnedWidth, binnedHeight, binnedChannels, binnedSlices, binnedFrames);
			makeRectangle(1, 1, binnedWidth - 2, binnedHeight - 2);
			getStatistics(area, mean, min, max, std);
			varianceOverBinningSq = (std * std) / (b * b);
			resultArrayX = Array.concat(resultArrayX, b);
			resultArrayY = Array.concat(resultArrayY, varianceOverBinningSq);
			if (varianceOverBinningSq > maxResultValue) {
				maxResultValue = varianceOverBinningSq;
			}
			close("DQE0_ActiveImage_Binned");
			Plot.create("DQE0 Noise Binning Results", "Binning", "Variance / Binning^2");
			Plot.setFrameSize(800, 600);
			Plot.setLimits(0, maxBinning, 0, maxResultValue * 1.05 );
			Plot.setColor("#1965B0");
			Plot.add("circles", resultArrayX, resultArrayY, "Data");
			Plot.update();
		}
		close("DQE0_ActiveImage");
		comparisons++;
	}
}

for (i = 1; i <= slices; i++) {
	close("DQE0_Image_" + i);
}

resultArray = newArray(0);
minBinning = Math.floor(maxBinning * 0.667);
for (i = 0; i < resultArrayX.length; i++) {
	if (resultArrayX[i] > minBinning) {
		resultArray = Array.concat(resultArray, resultArrayY[i]);
	}
}
Array.getStatistics(resultArray, min, max, mean, std);

meanX = newArray(-1, maxBinning + 1);
meanY = newArray(mean, mean);
stdError = std / Math.sqrt(comparisons);

dqe0 = (meanADU * meanADU) / (0.5 * mean * electronsPerPixel);
dqe0SigmaExposure = Math.sqrt(electronsPerPixel) * (meanADU * meanADU) / (0.5 * mean * electronsPerPixel * electronsPerPixel);
dqe0SigmaVariance = stdError * (meanADU * meanADU) / (0.5 * mean * mean * electronsPerPixel);
dqe0Error = Math.sqrt(dqe0SigmaExposure * dqe0SigmaExposure + dqe0SigmaVariance * dqe0SigmaVariance);

dqe0 = Math.round(dqe0 * 1000) / 1000.0;
dqe0Error = Math.round(dqe0Error * 1000) / 1000.0;

Plot.create("DQE0 Noise Binning Results", "Binning", "Variance / Binning^2");
Plot.setFrameSize(800, 600);
Plot.setLimits(0, maxBinning, 0, maxResultValue * 1.05 );
Plot.setColor("#CDD5EB", "#CDD5EB");
Plot.drawShapes("rectangles", -1, mean + stdError, maxBinning + 1, mean - stdError);
Plot.setColor("#1965B0");
Plot.add("circles", resultArrayX, resultArrayY, "Data");
Plot.add("line", meanX, meanY, "Mean");
Plot.addText("DQE(0) = " + dqe0 + " +/- " + dqe0Error, 0.02, 0.05);
Plot.update();

print("");
print("RESULTS");
print("Number of images : " + slices);
print("Mean intensity : " + meanADU);
print("ADU per electron : " + aduPerElectron);
print("Electrons per pixel : " + electronsPerPixel);
print("DQE(0) exposure sigma : " + dqe0SigmaExposure);
print("DQE(0) variance sigma : " + dqe0SigmaVariance);
print("DQE(0) : " + dqe0 + " +/- " + dqe0Error);

print("");
print("Done");
