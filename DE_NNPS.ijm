print("\\Clear");
print("DE_NNPS");
print("Version : 1.01");
print("Date : 2021.12.21");
print("Author : Benjamin Bammes, Direct Electron LP (bbammes@directelectron.com)");
print("License : GNU General Public License v2.0");
print("Requires : https://imagej.nih.gov/ij/plugins/radial-profile.html");
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

run("Duplicate...", "title=NNPS_CroppedStack duplicate");
Stack.getDimensions(width, height, channels, slices, frames);
if (width != height) {
	close("NNPS_CroppedStack");
	print("Error");
	exit("You must select a square ROI before running this macro.");
}

Stack.getStatistics(voxelCount, mean, min, max, std);
minValue = mean - std * 4.0;
maxValue = mean + std * 4.0;
run("Min...", "value=" + minValue + " stack");
run("Max...", "value=" + maxValue + " stack");

resultArrayX = newArray(0);
resultArrayY = newArray(0);

for (i = 1; i <= slices; i++) {
	Stack.setSlice(i);
	run("Set Label...", "label=NNPS_Image_" + i);
}
run("Stack to Images");

psWidth = 0;
comparisons = 0;

for (i = 1; i <= slices; i++) {
	for (j = i + 1; j <= slices; j++) {
		print("Processing frame " + i + " and " + j + "...");
		imageCalculator("Subtract create 32-bit", "NNPS_Image_" + i, "NNPS_Image_" + j);
		selectWindow("Result of NNPS_Image_" + i);
		rename("NNPS_ActiveImage");
		selectWindow("NNPS_ActiveImage");
		run("FFT Options...", "fft raw do");
		close("NNPS_ActiveImage");
		close("FFT of NNPS_ActiveImage");
		selectWindow("PS of NNPS_ActiveImage");
		getDimensions(psWidth, psHeight, psChannels, psSlices, psFrames);
		makeRectangle(0, 0, psWidth, psHeight);
		run("Radial Profile", "x=" + Math.floor(psWidth / 2) + " y=" + Math.floor(psHeight / 2) + " radius=" + Math.floor(psWidth / 2));
		Plot.getValues(plotX, plotY);
		for (n = 0; n < plotX.length; n++) {
			resultArrayX = Array.concat(resultArrayX, plotX[n] / (psWidth / 2.0));
			resultArrayY = Array.concat(resultArrayY, plotY[n]);
		}
		close("PS of NNPS_ActiveImage");
		close("Radial Profile Plot");
		comparisons++;
	}
}

for (i = 1; i <= slices; i++) {
	close("NNPS_Image_" + i);
}

print("Normalizing and averaging results...");

lowFreqArray = newArray(0);
for (i = 0; i < resultArrayX.length; i++) {
	if ((resultArrayX[i] >= 0.04) && (resultArrayX[i] <= 0.08)) {
		lowFreqArray = Array.concat(lowFreqArray, resultArrayY[i]);
	}
}
Array.getStatistics(lowFreqArray, lowFreqMin, lowFreqMax, lowFreqMean, lowFreqStd);
highFreqArray = newArray(0);
for (i = 0; i < resultArrayX.length; i++) {
	if ((resultArrayX[i] >= 0.80) && (resultArrayX[i] <= 0.96)) {
		highFreqArray = Array.concat(highFreqArray, resultArrayY[i]);
	}
}
Array.getStatistics(highFreqArray, highFreqMin, highFreqMax, highFreqMean, highFreqStd);
mean = lowFreqMean;
if (((lowFreqMean / highFreqMean) >= 0.666667) && ((lowFreqMean / highFreqMean) <= 1.5)) {
	print("  Counting mode detected");
	mean = highFreqMean;
}
for (i = 0; i < resultArrayX.length; i++) {
	resultArrayY[i] = resultArrayY[i] / mean;
}

print("NPS(0) : " + mean);

nnpsLength = Math.floor(psWidth / 2);

xArray = newArray(0);
meanArray = newArray(0);
minArray = newArray(0);
maxArray = newArray(0);
for (n = 0; n < nnpsLength; n++) {
	valuesArray = newArray(0);
	for (i = 0; i < resultArrayX.length; i++) {
		if ((resultArrayX[i] >= ((n - 0.5) / nnpsLength)) && (resultArrayX[i] < ((n + 0.5) / nnpsLength))) {
			valuesArray = Array.concat(valuesArray, resultArrayY[i]);
		}
	}
	if (valuesArray.length > 0) {
		Array.getStatistics(valuesArray, valuesMin, valuesMax, valuesMean, valuesStd);
		valuesStd = valuesStd / comparisons;
		xArray = Array.concat(xArray, (n + 0.0) / nnpsLength);
		meanArray = Array.concat(meanArray, valuesMean);
		minArray = Array.concat(minArray, valuesMean - valuesStd / Math.sqrt(comparisons));
		maxArray = Array.concat(maxArray, valuesMean + valuesStd / Math.sqrt(comparisons));
	}
}

smoothingWidth = Math.floor(nnpsLength / 16);
smoothingHalfWidth = Math.floor(smoothingWidth / 2);

equationToFit = "y = a*exp(b*x)";

meanArray[0] = 1.0;
minArray[0] = 1.0;
maxArray[0] = 1.0;

fitArrayX = newArray(0);
fitArrayMeanY = newArray(0);
fitArrayMinY = newArray(0);
fitArrayMaxY = newArray(0);
for (n = 0; n < smoothingHalfWidth; n++) {
	if ((meanArray[n] > 0.8) && (meanArray[n] < 1.2)) {
		fitArrayX = Array.concat(fitArrayX, n);
		fitArrayMeanY = Array.concat(fitArrayMeanY, meanArray[n]);
		fitArrayMinY = Array.concat(fitArrayMinY, minArray[n]);
		fitArrayMaxY = Array.concat(fitArrayMaxY, maxArray[n]);
	}	
}
Fit.doFit(equationToFit, fitArrayX, fitArrayMeanY, newArray(meanArray[0], -0.01));
for (n = 0; n < smoothingHalfWidth; n++) {
	meanArray[n] = Fit.f(n + 0.0);
}
Fit.doFit(equationToFit, fitArrayX, fitArrayMinY, newArray(minArray[0], -0.01));
for (n = 0; n < smoothingHalfWidth; n++) {
	minArray[n] = Fit.f(n + 0.0);
}
Fit.doFit(equationToFit, fitArrayX, fitArrayMaxY, newArray(maxArray[0], -0.01));
for (n = 0; n < smoothingHalfWidth; n++) {
	maxArray[n] = Fit.f(n + 0.0);
}

for (n = smoothingHalfWidth; n < meanArray.length - smoothingHalfWidth; n++) {
	fitArrayX = newArray(0);
	fitArrayMeanY = newArray(0);
	fitArrayMinY = newArray(0);
	fitArrayMaxY = newArray(0);
	for (k = (-1 * smoothingHalfWidth); k <= smoothingHalfWidth; k++) {
		fitArrayX = Array.concat(fitArrayX, k);
		fitArrayMeanY = Array.concat(fitArrayMeanY, meanArray[n + k]);
		fitArrayMinY = Array.concat(fitArrayMinY, minArray[n + k]);
		fitArrayMaxY = Array.concat(fitArrayMaxY, maxArray[n + k]);
	}
	Fit.doFit(equationToFit, fitArrayX, fitArrayMeanY, newArray(meanArray[n], -0.5));
	meanArray[n] = Fit.f(0.0);
	Fit.doFit(equationToFit, fitArrayX, fitArrayMinY, newArray(minArray[n], -0.5));
	minArray[n] = Fit.f(0.0);
	Fit.doFit(equationToFit, fitArrayX, fitArrayMaxY, newArray(maxArray[n], -0.5));
	maxArray[n] = Fit.f(0.0);
}

fitArrayX = newArray(0);
fitArrayMeanY = newArray(0);
fitArrayMinY = newArray(0);
fitArrayMaxY = newArray(0);
for (n = meanArray.length - smoothingHalfWidth; n < meanArray.length; n++) {
	fitArrayX = Array.concat(fitArrayX, n - (meanArray.length - smoothingHalfWidth));
	fitArrayMeanY = Array.concat(fitArrayMeanY, meanArray[n]);
	fitArrayMinY = Array.concat(fitArrayMinY, minArray[n]);
	fitArrayMaxY = Array.concat(fitArrayMaxY, maxArray[n]);
}
Fit.doFit(equationToFit, fitArrayX, fitArrayMeanY, newArray(meanArray[meanArray.length - smoothingHalfWidth], -0.01));
for (n = meanArray.length - smoothingHalfWidth; n < meanArray.length; n++) {
	meanArray[n] = Fit.f(n - (meanArray.length - smoothingHalfWidth));
}
Fit.doFit(equationToFit, fitArrayX, fitArrayMinY, newArray(minArray[meanArray.length - smoothingHalfWidth], -0.01));
for (n = meanArray.length - smoothingHalfWidth; n < meanArray.length; n++) {
	minArray[n] = Fit.f(n - (meanArray.length - smoothingHalfWidth));
}
Fit.doFit(equationToFit, fitArrayX, fitArrayMaxY, newArray(maxArray[meanArray.length - smoothingHalfWidth], -0.01));
for (n = meanArray.length - smoothingHalfWidth; n < meanArray.length; n++) {
	maxArray[n] = Fit.f(n - (meanArray.length - smoothingHalfWidth));
}

for (n = smoothingHalfWidth - 2; n < smoothingHalfWidth + 1; n++) {
	meanArray[n] = (meanArray[n - 2] + meanArray[n - 1] + meanArray[n + 1] + meanArray[n + 2]) / 4.0;
	minArray[n] = (minArray[n - 2] + minArray[n - 1] + minArray[n + 1] + minArray[n + 2]) / 4.0;
	maxArray[n] = (maxArray[n - 2] + maxArray[n - 1] + maxArray[n + 1] + maxArray[n + 2]) / 4.0;
}

Plot.create("NNPS Results", "1 / Nyquist", "NNPS");
Plot.setFrameSize(800, 600);
Plot.setLimits(0, 1, 0, 1.1 );
Plot.setColor("#CDD5EB");
Plot.add("line", xArray, minArray, "Mean - StdErr");
Plot.add("line", xArray, maxArray, "Mean + StdErr");
Plot.setColor("#1965B0");
Plot.add("circles", resultArrayX, resultArrayY, "Data");
Plot.add("line", xArray, meanArray, "Mean");
Plot.show();

print("");
print("Done");
