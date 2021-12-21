print("\\Clear");
print("DE_MagCal");
print("Version : 1.01");
print("Date : 2021.12.21");
print("Author : Benjamin Bammes, Direct Electron LP (bbammes@directelectron.com)");
print("License : GNU General Public License v2.0");
print("Requires : https://imagej.nih.gov/ij/plugins/radial-profile.html");
print("");

openImages = getList("window.titles");
if (openImages.length < 1) {
	print("Error");
	exit("No image found.");
}
if (getValue("selection.size") != 4) {
	print("Error");
	exit("You must select a single square ROI before running this macro.");
}
activeImageTitle = getTitle();

apix = 0.0;

typeOfData = getNumber("Type of Data (0 = crossed-lines (low magnifications),  1 = gold (high magnifications))", 0);

if (typeOfData == 0) {
	print("Finding pixel size using crossed-lines (for low magnifications)");
}
else {
	print("Finding pixel size using gold lattice spacing (for high magnifications)");
}

print("Calculating FFT...");

run("FFT Options...", "fft raw do");
close("FFT of " + activeImageTitle);
selectWindow("PS of " + activeImageTitle);
getDimensions(psWidth, psHeight, psChannels, psSlices, psFrames);
if ((psWidth < 512) || (psHeight < 512)) {
	close("PS of " + activeImageTitle);
	print("Error");
	exit("ROI size must be at least 512 x 512.");
}
rename("MagCal_PS_Full");

if (typeOfData == 0) {

	setPixel(psWidth / 2, psHeight / 2, 0.0);
	
	makeRectangle(psWidth * 3 / 8, psHeight * 3 / 8, psWidth / 4, psHeight / 4);
	run("Crop");
	rename("MagCal_PS");
	getDimensions(psWidthCropped, psHeightCropped, psChannelsCropped, psSlicesCropped, psFramesCropped);
	getStatistics(psArea, psMean, psMin, psMax, psStd, psHistogram);

	psBGMaxValue = 0.0;
	for (dx = -4; dx <= 4; dx++) {
		for (dy = -4; dy <= 4; dy++) {
			if ((dx == 0) && (dy == 0)) {
				continue;
			}
			psValue = getPixel(dx + psWidthCropped / 2, dy + psHeightCropped / 2);
			if (psValue > psBGMaxValue) {
				psBGMaxValue = psValue;
			}
		}
	}

	print("Finding FFT peaks...");

	psPeaksMaxNumber = 1000;
	psPeaksX = newArray(psPeaksMaxNumber);
	psPeaksY = newArray(psPeaksMaxNumber);
	psPeaksNumber = 0;
	for (x = 1; x < (psWidthCropped - 1); x++) {
		if (psPeaksNumber >= psPeaksMaxNumber) {
			break;
		}
		for (y = 1; y < (psHeightCropped - 1); y++) {
			if (psPeaksNumber >= psPeaksMaxNumber) {
				break;
			}
			psValue = getPixel(x, y);
			if (psValue > psBGMaxValue) {
				if (psValue > getPixel(x - 1, y)) {
					if (psValue > getPixel(x + 1, y)) {
						if (psValue > getPixel(x, y - 1)) {
							if (psValue > getPixel(x, y + 1)) {
								if (psValue > getPixel(x - 1, y - 1)) {
									if (psValue > getPixel(x - 1, y + 1)) {
										if (psValue > getPixel(x + 1, y - 1)) {
											if (psValue > getPixel(x + 1, y + 1)) {
												
													x1 = x - 1.0;
													x2 = x + 0.0;
													x3 = x + 1.0;
													y1 = getPixel(x - 1, y);
													y2 = getPixel(x, y);
													y3 = getPixel(x + 1, y);
													denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
													fitA = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
													fitB = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
													fitC = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
													bestX = -1.0 * fitB / ( 2 * fitA);
												
													x1 = y - 1.0;
													x2 = y + 0.0;
													x3 = y + 1.0;
													y1 = getPixel(x, y - 1);
													y2 = getPixel(x, y);
													y3 = getPixel(x, y + 1);
													denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
													fitA = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
													fitB = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
													fitC = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
													bestY = -1.0 * fitB / ( 2 * fitA);

													psPeaksX[psPeaksNumber] = bestX;
													psPeaksY[psPeaksNumber] = bestY;
													psPeaksNumber++;
													
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}

	if (psPeaksNumber < 1) {
		print("Error");
		exit("Could not find peaks in the FFT.");
	}

	print("  Found " + psPeaksNumber + " FFT peaks");

	print("Calculating distances between FFT peaks...");

	minValidDistance = 6.0;

	distancesFromEachOther = newArray(psPeaksNumber * psPeaksNumber);
	distancesFromCenter = newArray(psPeaksNumber);
	distancesNumber = 0;
	for (i = 0; i < psPeaksNumber; i++) {
		distancesFromCenter[i] = Math.sqrt((psPeaksX[i] - psWidthCropped / 2.0) * (psPeaksX[i] - psWidthCropped / 2.0) + (psPeaksY[i] - psHeightCropped / 2.0) * (psPeaksY[i] - psHeightCropped / 2.0));
		if (distancesFromCenter[i] < minValidDistance) {
			distancesFromCenter[i] = 0.0;
			continue;
		}
		for (j = i + 1; j < psPeaksNumber; j++) {
			distancesFromEachOther[distancesNumber] = Math.sqrt((psPeaksX[i] - psPeaksX[j]) * (psPeaksX[i] - psPeaksX[j]) + (psPeaksY[i] - psPeaksY[j]) * (psPeaksY[i] - psPeaksY[j]));
			distancesNumber++;
		}
	}

	for (i = 0; i < distancesNumber; i++) {
		inDistancesFromCenterArray = false;
		for (j = 0; j < psPeaksNumber; j++) {
			if ((distancesFromCenter[j] >= minValidDistance) && (Math.abs(distancesFromEachOther[i] - distancesFromCenter[j]) < minValidDistance)) {
				inDistancesFromCenterArray = true;
				break;
			}
		}
		if (!inDistancesFromCenterArray) {
			distancesFromEachOther[i] = 0.0;
		}
	}

	print("Finding most common minimum distance...");
	
	distancesClassifiedValues = newArray(distancesNumber);
	distancesClassifiedCounts = newArray(distancesNumber);
	distancesClassifiedNumber = 0;
	for (i = 0; i < distancesNumber; i++) {
		if (distancesFromEachOther[i] < minValidDistance) {
			continue;
		}
		notClassified = true;
		for (j = 0; j < distancesClassifiedNumber; j++) {
			if (Math.abs(distancesFromEachOther[i] - distancesClassifiedValues[j]) < minValidDistance) {
				distancesClassifiedValues[j] *= distancesClassifiedCounts[j];
				distancesClassifiedValues[j] += distancesFromEachOther[i];
				distancesClassifiedCounts[j] += 1.0;
				distancesClassifiedValues[j] /= distancesClassifiedCounts[j];
				notClassified = false;
				break;
			}
		}
		if (notClassified) {
			distancesClassifiedValues[distancesClassifiedNumber] += distancesFromEachOther[i];
			distancesClassifiedCounts[distancesClassifiedNumber] += 1.0;
			distancesClassifiedNumber++;
		}
	}

	minDistance = psWidth;
	for (j = 0; j < distancesClassifiedNumber; j++) {
		if ((distancesClassifiedValues[j] < minDistance) && (distancesClassifiedCounts[j] > 1.0)) {
			minDistance = distancesClassifiedValues[j];
		}
	}

	print("  Distances from each other : ");
	for (j = 0; j < distancesClassifiedNumber; j++) {
		print("    " + distancesClassifiedValues[j]);
	}

	print("  Fourier-space distance : " + minDistance + " pixels");
	
	fractionOfNyquist = minDistance / (psWidth / 2.0);
	lineSpacingInAngstroms = 10000000.0 / 2160.0;
	apix = round(10000.0 * lineSpacingInAngstroms * fractionOfNyquist) / 10000.0;

	makeRectangle((psWidthCropped / 2.0 - minDistance * 4.0), Math.round(psHeightCropped / 2.0 - minDistance * 4.0), Math.round(8.0 * minDistance), Math.round(8.0 * minDistance));
	run("Crop");
	getDimensions(psWidthCropped, psHeightCropped, psChannelsCropped, psSlicesCropped, psFramesCropped);
	run("Min...", "value=" + psBGMaxValue);
	run("Max...", "value=" + psBGMaxValue * 2.0);
	run("Enhance Contrast", "saturated=0.75");

	makeOval(Math.round(psWidthCropped / 2.0 - minDistance), Math.round(psHeightCropped / 2.0 - minDistance), Math.round(2.0 * minDistance), Math.round(2.0 * minDistance));
	
}

else if (typeOfData == 1) {

	run("Scale...", "x=.25 y=.25 width=" + (psWidth / 4) + " height=" + (psHeight / 4) + " interpolation=Bilinear average create");
	rename("MagCal_PS");
	close("MagCal_PS_Full");
	
	selectWindow("MagCal_PS");
	getDimensions(psWidth, psHeight, psChannels, psSlices, psFrames);
	
	makeRectangle(psWidth * 17 / 32, psHeight * 17 / 32, psWidth / 3, psHeight / 3);
	run("Enhance Contrast", "saturated=0.75");
	
	run("Select All");
	run("Radial Profile", "x=" + Math.floor(psWidth / 2) + " y=" + Math.floor(psHeight / 2) + " radius=" + Math.floor(psWidth / 2));
	Plot.getValues(plotX, plotY);
	close("Radial Profile Plot");
	
	print("Finding gold peak...");
	
	bgSubtractedY = newArray(plotX.length);
	for (n = 0; n < plotX.length; n++) {
		if ((n < maxOf(8, Math.round(plotX.length / 16.0))) || (n >= (plotX.length - 8))) {
			bgSubtractedY[n] = 0.0;
		}
		else {
			bg = 0.0;
			bgN = 0.0;
			for (i = -8; i <= 8; i++) {
				if (i != 0) {
					bg += plotY[n + i];
					bgN += 1.0;
				}
			}
			bgSubtractedY[n] = plotY[n] - (bg / bgN);
		}
	}
	
	maxValue = 0.0;
	maxLocation = 0;
	for (n = 0; n < plotX.length; n++) {
		if (bgSubtractedY[n] > maxValue) {
			maxValue = bgSubtractedY[n];
			maxLocation = n;
		}
	}
	
	x1 = plotX[maxLocation - 1];
	x2 = plotX[maxLocation];
	x3 = plotX[maxLocation + 1];
	y1 = bgSubtractedY[maxLocation - 1];
	y2 = bgSubtractedY[maxLocation];
	y3 = bgSubtractedY[maxLocation + 1];
	denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
	fitA = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
	fitB = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
	fitC = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
	peakLocation = -1.0 * fitB / ( 2 * fitA);
	
	peakX = newArray(1);
	peakX[0] = peakLocation;
	peakY = newArray(1);
	peakY[0] = maxValue;
	
	print("  Peak found at radius " + peakX[0]);

	makeOval(Math.round(psWidth / 2.0 - peakLocation), Math.round(psHeight / 2.0 - peakLocation), Math.round(2.0 * peakLocation), Math.round(2.0 * peakLocation));
	
	apix = round(10000.0 * (peakLocation / (psWidth / 2.0)) * 2.321 / 2.0) / 10000.0;

}

else {
	print("Error");
	exit("Invalid data type.");
}

print("");
print("RESULTS");
print("Angstroms/pixel : " + apix);

print("");
print("Done");
