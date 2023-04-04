print("\\Clear");
print("DE_FalsePos");
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

thresholdArray = newArray(0);
epsArray = newArray(0);
validityArray = newArray(0);

print("Calculating the false positive rate...");

Stack.getDimensions(width, height, channels, slices, frames);
for (i = 1; i <= slices; i++) {
	Stack.setSlice(i);
	label = getInfo("slice.label");
	indexOfFPSubstring = label.indexOf("_FP_");
	if (indexOfFPSubstring < 1) {
		print("Error");
		exit("The label of each image must include '_FP_' followed by a three digit number representing the sensor threshold.");
	}
	thresholdArray = Array.concat(thresholdArray, parseInt(label.substring(indexOfFPSubstring + 4, indexOfFPSubstring + 7)));
	getStatistics(selectionArea, selectionMean);
	epsArray = Array.concat(epsArray, selectionMean / (aduPerE * exposureTime));
	validityArray = Array.concat(validityArray, 1);
}

print("Sorting arrays...");
thresholdArraySorted = newArray(0);
epsArraySorted = newArray(0);
elementsRemaining = epsArray.length;
while (elementsRemaining > 0) {
	nextThreshold = 0;
	nextIndex = -1;
	for (i = 0; i < epsArray.length; i++) {
		if (validityArray[i] > 0) {
			if (nextIndex < 0) {
				nextThreshold = thresholdArray[i];
				nextIndex = i;
			}
			else if (thresholdArray[i] < nextThreshold) {
				nextThreshold = thresholdArray[i];
				nextIndex = i;
			}
		}
	}
	if (nextIndex >= 0) {
		thresholdArraySorted = Array.concat(thresholdArraySorted, thresholdArray[nextIndex]);
		epsArraySorted = Array.concat(epsArraySorted, epsArray[nextIndex]);
		validityArray[nextIndex] = 0;
		elementsRemaining -= 1;
	}
	else {
		break;
	}
}

print("Generating plots...");
Array.getStatistics(thresholdArraySorted, thresholdMin, thresholdMax);
Array.getStatistics(epsArraySorted, epsMin, epsMax);
Plot.create("False Positive Rate Results", "Sensor Threshold", "False Positive Rate (eps)");
Plot.setFrameSize(800, 600);
Plot.setLimits(thresholdMin, thresholdMax, 0, epsMax * 1.1 );
Plot.setColor("#1965B0");
Plot.add("circles", thresholdArraySorted, epsArraySorted, "Data");
Plot.add("line", thresholdArraySorted, epsArraySorted, "Line");
Plot.show();

thresholdRecommended = 99999;
for (i = 0; i < epsArraySorted.length; i++) {
	if (epsArraySorted[i] < 0.01) {
		if (thresholdRecommended > thresholdArraySorted[i]) {
			thresholdRecommended = thresholdArraySorted[i];
		}
	}
}

print("");
print("RESULTS");
if (thresholdRecommended < 99990)
	print("Recommended Threshold : " + thresholdRecommended);
else
	print("Recommended Threshold : Unknown");

print("");
print("Done");
