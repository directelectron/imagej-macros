print("\\Clear");
print("DE_CoincLoss");
print("Version : 1.00");
print("Date : 2023.03.30");
print("Author : Benjamin Bammes, Direct Electron LP (bbammes@directelectron.com)");
print("License : GNU General Public License v2.0");
print("");

openImages = getList("image.titles");
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

aduPerE = getNumber("ADU per primary electron", 0);
if (aduPerE == 0) {
	print("Error");
	exit("The value of a single electron cannot be zero.");
}

exposureTime = getNumber("Exposure time of each image (seconds)", 1);
if (exposureTime <= 0) {
	print("Error");
	exit("The value of exposure time must be positive.");
}

pixelSizeArray = newArray(0);
epsArray = newArray(0);
validityArray = newArray(0);

print("Calculating the intensity of each image...");

Stack.getDimensions(width, height, channels, slices, frames);
for (i = 1; i <= slices; i++) {
	Stack.setSlice(i);
	label = getInfo("slice.label");
	indexOfCLSubstring = label.indexOf("_CL_");
	if (indexOfCLSubstring < 1) {
		print("Error");
		exit("The label of each image must include '_CL_' followed by a four digit number representing the pixel size in picometers.");
	}
	pixelSizeArray = Array.concat(pixelSizeArray, parseInt(label.substring(indexOfCLSubstring + 4, indexOfCLSubstring + 8)));
	getStatistics(selectionArea, selectionMean);
	epsArray = Array.concat(epsArray, selectionMean / (aduPerE * exposureTime));
	validityArray = Array.concat(validityArray, 1);
}

print("Sorting arrays...");
pixelSizeArraySorted = newArray(0);
epsArraySorted = newArray(0);
elementsRemaining = epsArray.length;
while (elementsRemaining > 0) {
	nextPixelSize = 0;
	nextIndex = -1;
	for (i = 0; i < epsArray.length; i++) {
		if (validityArray[i] > 0) {
			if (nextIndex < 0) {
				nextPixelSize = pixelSizeArray[i];
				nextIndex = i;
			}
			else if (pixelSizeArray[i] < nextPixelSize) {
				nextPixelSize = pixelSizeArray[i];
				nextIndex = i;
			}
		}
	}
	if (nextIndex >= 0) {
		pixelSizeArraySorted = Array.concat(pixelSizeArraySorted, pixelSizeArray[nextIndex]);
		epsArraySorted = Array.concat(epsArraySorted, epsArray[nextIndex]);
		validityArray[nextIndex] = 0;
		elementsRemaining -= 1;
	}
	else {
		break;
	}
}

print("Calculating coincidence loss...");
inputEPSArraySorted = newArray(0);
for (i = 0; i < epsArraySorted.length; i++) {
	inputEPSArraySorted = Array.concat(inputEPSArraySorted, epsArraySorted[0] * pixelSizeArraySorted[i] * pixelSizeArraySorted[i] / (pixelSizeArraySorted[0] * pixelSizeArraySorted[0]));
}

Array.getStatistics(inputEPSArraySorted, inputEPSMin, inputEPSMax);
Array.getStatistics(epsArraySorted, epsMin, epsMax);
maxEPSToFit = Math.ceil(inputEPSMax / 5) * 5;

equationToFit = "y = (1-exp(-1*a*x))/a";
initialGuesses = newArray(1);
initialGuesses[0] = 0.01;
Fit.doFit(equationToFit, inputEPSArraySorted, epsArraySorted, initialGuesses);
curveInputEPS = newArray(0);
curveOutputEPS = newArray(0);
for (eps = 0; eps < maxEPSToFit; eps++) {
	curveInputEPS = Array.concat(curveInputEPS, eps);
	curveOutputEPS = Array.concat(curveOutputEPS, Fit.f(eps));
}

eps = 0.1;
epsAt05CL = -1.0;
epsAt10CL = -1.0;
epsAt20CL = -1.0;
while (epsAt20CL < 0.0) {
	clFraction = 1.0 - (Fit.f(eps) / eps);
	if (epsAt05CL < 0.0) {
		if (clFraction > 0.05) {
			epsAt05CL = eps;
		}
	}
	if (epsAt10CL < 0.0) {
		if (clFraction > 0.10) {
			epsAt10CL = eps;
		}
	}
	if (epsAt20CL < 0.0) {
		if (clFraction > 0.20) {
			epsAt20CL = eps;
		}
	}
	eps += 0.1;
	if (eps > 1000.0) {
		break;
	}
}

print("Generating plots...");
Plot.create("False Positive Rate Results", "Sensor Threshold", "False Positive Rate (eps)");
Plot.setFrameSize(800, 600);
Plot.setLimits(0, maxEPSToFit, 0, maxEPSToFit );
Plot.setColor("#1965B0");
Plot.add("circles", inputEPSArraySorted, epsArraySorted, "Data");
Plot.add("line", curveInputEPS, curveOutputEPS, "Line");
Plot.show();

print("");
print("RESULTS");
print("y = (1-exp(-1*a*x))/a, where a = " + Fit.p(0));
print("Rsquared = " + Fit.rSquared);
if (Fit.rSquared < 0.995)
	print("The quality of fit is less than expected.\nRemove any images that appear to be obvious outliers (especially at high exposure rate) and try again.");
if (epsAt05CL > 0.0)
	print("5% coincidence loss at " + epsAt05CL + " eps");
else
	print("5% coincidence loss at unknown");
if (epsAt10CL > 0.0)
	print("10% coincidence loss at " + epsAt10CL + " eps");
else
	print("10% coincidence loss at unknown");
if (epsAt20CL > 0.0)
	print("20% coincidence loss at " + epsAt20CL + " eps");
else
	print("20% coincidence loss at unknown");

print("");
print("Done");
