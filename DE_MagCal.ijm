print("\\Clear");
print("DE_MagCal");
print("Version : 1.03");
print("Date : 2022.08.16");
print("Author : Benjamin Bammes, Direct Electron LP (bbammes@directelectron.com)");
print("License : GNU General Public License v2.0");
print("Requires : https://imagej.nih.gov/ij/plugins/radial-profile.html");
print("");

openImages = getList("image.titles");
if (openImages.length < 1) {
	print("Error");
	exit("No image found.");
}
if (getValue("selection.size") != 4) {
	run("Select All");
}
activeImageTitle = getTitle();

apix = 0.0;

typeOfData = getNumber("Type of Data (0 = crossed-lines (low magnifications),  1 = gold (high magnifications))", 0);

useEstimatedAngle = false;
estimatedAngle = 0.0;

if (typeOfData == 0) {
	print("Finding pixel size using crossed-lines (for low magnifications)");
	useEstimatedAngle = getBoolean("Do you want to enter an estimated angle for the line grating?");
	if (useEstimatedAngle) {
		estimatedAngle = getNumber("Estimated angle in degrees of the line grating", 0);
		if (estimatedAngle < 180.0) {
			estimatedAngle += 180.0;
		}
		if (estimatedAngle >= 90.0) {
			estimatedAngle -= 90.0;
		}
		estimatedAngle *= PI / 180.0;
	}
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

	for (x = 0; x < psWidth; x++) {
		setPixel(x, psHeight / 2, 0.0);
	}
	for (y = 0; y < psHeight; y++) {
		setPixel(psWidth / 2, y, 0.0);
	}
	
	makeRectangle(psWidth * 7 / 16, psHeight * 7 / 16, psWidth / 8, psHeight / 8);
	run("Crop");
	rename("MagCal_PS_Raw");
	getDimensions(psWidthCropped, psHeightCropped, psChannelsCropped, psSlicesCropped, psFramesCropped);

	run("Duplicate...", " ");
	rename("MagCal_PS_BG");
	run("Median...", "radius=8");
	imageCalculator("Subtract create 32-bit", "MagCal_PS_Raw", "MagCal_PS_BG");
	close("MagCal_PS_Raw");
	close("MagCal_PS_BG");
	selectWindow("Result of MagCal_PS_Raw");
	rename("MagCal_PS");
	getStatistics(psArea, psMean, psMin, psMax, psStd, psHistogram);
	psBGMaxValue = psMean * 10.0;
	run("Enhance Contrast", "saturated=0.75");

	print("Finding FFT peaks...");

	psPeaksMaxNumber = 500;
	psPeaksX = newArray(psPeaksMaxNumber);
	psPeaksY = newArray(psPeaksMaxNumber);
	psPeaksNumber = 0;
	for (x = 4; x < (psWidthCropped - 4); x++) {
		if (psPeaksNumber >= psPeaksMaxNumber) {
			break;
		}
		for (y = 4; y < (psHeightCropped - 4); y++) {
			if (psPeaksNumber >= psPeaksMaxNumber) {
				break;
			}
			distanceFromCenter = Math.sqrt((x - psWidthCropped / 2.0) * (x - psWidthCropped / 2.0) + (y - psHeightCropped / 2.0) * (y - psHeightCropped / 2.0));
			if (distanceFromCenter <= 16) {
				continue;
			}
			psValue = getPixel(x, y);
			if (psValue > psBGMaxValue) {
				
				isLocalMax = true;
				for (dx = -4; dx <= 4; dx++) {
					for (dy = -4; dy <= 4; dy++) {
						if (Math.sqrt(dx * dx + dy * dy) <= 4) {
							if (psValue < getPixel(x + dx, y + dy)) {
								isLocalMax = false;
								break;
							}
						}
					}
					if (!isLocalMax) {
						break;
					}
				}
				
				if (isLocalMax) {
													
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

	if (true) {
		run("Duplicate...", " ");
		rename("MagCal_PS_Peaks");
		run("Multiply...", "value=0.0");
		for (i = 0; i < psPeaksNumber; i++) {
			setPixel(psPeaksX[i], psPeaksY[i], 1.0);
		}
		selectWindow("MagCal_PS");
	}

	if (psPeaksNumber < 1) {
		print("Error");
		exit("Could not find peaks in the FFT.");
	}

	print("  Found " + psPeaksNumber + " FFT peaks");

	print("Calculating distances between FFT peaks...");

	minValidAngularDifference = PI / 16;
	minValidDistance = 2.0;

	anglesFromEachOther = newArray(psPeaksNumber * psPeaksNumber);
	distancesFromEachOther = newArray(psPeaksNumber * psPeaksNumber);
	distancesNumber = 0;
	for (i = 0; i < psPeaksNumber; i++) {
		for (j = i + 1; j < psPeaksNumber; j++) {
			angleFromEachOther = Math.atan2(psPeaksX[i] - psPeaksX[j], psPeaksY[i] - psPeaksY[j]);
			if (angleFromEachOther < 0.0) {
				angleFromEachOther += PI;
			}
			if (angleFromEachOther >= (PI / 2.0)) {
				angleFromEachOther -= PI / 2.0;
			}
			if (useEstimatedAngle) {
				if ((Math.abs(estimatedAngle - angleFromEachOther) > minValidAngularDifference) && (Math.abs(estimatedAngle - angleFromEachOther - PI / 2.0) > minValidAngularDifference) && (Math.abs(estimatedAngle - angleFromEachOther + PI / 2.0) > minValidAngularDifference)) {
					continue;
				}
			}
			distanceFromEachOther = Math.sqrt((psPeaksX[i] - psPeaksX[j]) * (psPeaksX[i] - psPeaksX[j]) + (psPeaksY[i] - psPeaksY[j]) * (psPeaksY[i] - psPeaksY[j]));
			if (distanceFromEachOther > minValidDistance) {
				anglesFromEachOther[distancesNumber] = angleFromEachOther;
				distancesFromEachOther[distancesNumber] = distanceFromEachOther;
				distancesNumber++;
			}
		}
	}

	print("Grouping and averaging distances...");
	
	anglesClassifiedValues = newArray(distancesNumber);
	distancesClassifiedValues = newArray(distancesNumber);
	distancesClassifiedCounts = newArray(distancesNumber);
	distancesClassifiedNumber = 0;
	for (i = 0; i < distancesNumber; i++) {
		notClassified = true;
		for (j = 0; j < distancesClassifiedNumber; j++) {
			if ((Math.abs(anglesFromEachOther[i] - anglesClassifiedValues[j]) < minValidAngularDifference) || (Math.abs(anglesFromEachOther[i] - anglesClassifiedValues[j] - PI / 2.0) < minValidAngularDifference) || (Math.abs(anglesFromEachOther[i] - anglesClassifiedValues[j] + PI / 2.0) < minValidAngularDifference)) {
				if (Math.abs(distancesFromEachOther[i] - distancesClassifiedValues[j]) < minValidDistance) {
					anglesClassifiedValues[j] *= distancesClassifiedCounts[j];
					anglesClassifiedValues[j] += anglesFromEachOther[i];
					distancesClassifiedValues[j] *= distancesClassifiedCounts[j];
					distancesClassifiedValues[j] += distancesFromEachOther[i];
					distancesClassifiedCounts[j] += 1.0;
					anglesClassifiedValues[j] /= distancesClassifiedCounts[j];
					distancesClassifiedValues[j] /= distancesClassifiedCounts[j];
					notClassified = false;
					break;
				}
			}
		}
		if (notClassified) {
			anglesClassifiedValues[distancesClassifiedNumber] += anglesFromEachOther[i];
			distancesClassifiedValues[distancesClassifiedNumber] += distancesFromEachOther[i];
			distancesClassifiedCounts[distancesClassifiedNumber] += 1.0;
			distancesClassifiedNumber++;
		}
	}

	bestAngles = newArray(0, 0, 0, 0);
	bestDistances = newArray(4096.0, 4096.0, 4096.0, 4096.0);
	bestDistancesCount = newArray(0, 0, 0, 0);
	for (j = 0; j < distancesClassifiedNumber; j++) {
		if ((distancesClassifiedCounts[j] > 64) || (distancesClassifiedCounts[j] > bestDistancesCount[0]) || (distancesClassifiedCounts[j] > bestDistancesCount[1]) || (distancesClassifiedCounts[j] > bestDistancesCount[2]) || (distancesClassifiedCounts[j] > bestDistancesCount[3])) {
			if (distancesClassifiedValues[j] < bestDistances[0]) {
				bestAngles[3] = bestAngles[2];
				bestDistances[3] = bestDistances[2];
				bestDistancesCount[3] = bestDistancesCount[2];
				bestAngles[2] = bestAngles[1];
				bestDistances[2] = bestDistances[1];
				bestDistancesCount[2] = bestDistancesCount[1];
				bestAngles[1] = bestAngles[0];
				bestDistances[1] = bestDistances[0];
				bestDistancesCount[1] = bestDistancesCount[0];
				bestAngles[0] = anglesClassifiedValues[j] * 180.0 / PI;
				bestDistances[0] = distancesClassifiedValues[j];
				bestDistancesCount[0] = distancesClassifiedCounts[j];
			}
			else if (distancesClassifiedValues[j] < bestDistances[1]) {
				bestAngles[3] = bestAngles[2];
				bestDistances[3] = bestDistances[2];
				bestDistancesCount[3] = bestDistancesCount[2];
				bestAngles[2] = bestAngles[1];
				bestDistances[2] = bestDistances[1];
				bestDistancesCount[2] = bestDistancesCount[1];
				bestAngles[1] = anglesClassifiedValues[j] * 180.0 / PI;
				bestDistances[1] = distancesClassifiedValues[j];
				bestDistancesCount[1] = distancesClassifiedCounts[j];
			}
			else if (distancesClassifiedValues[j] < bestDistances[2]) {
				bestAngles[3] = bestAngles[2];
				bestDistances[3] = bestDistances[2];
				bestDistancesCount[3] = bestDistancesCount[2];
				bestAngles[2] = anglesClassifiedValues[j] * 180.0 / PI;
				bestDistances[2] = distancesClassifiedValues[j];
				bestDistancesCount[2] = distancesClassifiedCounts[j];
			}
			else if (distancesClassifiedValues[j] < bestDistances[3]) {
				bestAngles[3] = anglesClassifiedValues[j] * 180.0 / PI;
				bestDistances[3] = distancesClassifiedValues[j];
				bestDistancesCount[3] = distancesClassifiedCounts[j];
			}
		}
	}

	print("  Distances from each other : ");
	bestDistance = bestDistances[0];
	for (i = 0; i < 4; i++) {
		if (bestDistances[i] == bestDistance) {
			print("    " + bestDistances[i] + " pixels, " + bestAngles[i] + " deg (" + bestDistancesCount[i] + ") - BEST");
		} else {
			print("    " + bestDistances[i] + " pixels, " + bestAngles[i] + " deg (" + bestDistancesCount[i] + ")");
		}
	}
	
	fractionOfNyquist = bestDistance / (psWidth / 2.0);
	lineSpacingInAngstroms = 10000000.0 / 2160.0;
	apix = round(10000.0 * lineSpacingInAngstroms * fractionOfNyquist * 0.5) / 10000.0;

	selectWindow("MagCal_PS");
	run("Min...", "value=" + psBGMaxValue);
	run("Max...", "value=" + psBGMaxValue * 10.0);
	run("Enhance Contrast", "saturated=0.75");

	makeOval(Math.round(psWidthCropped / 2.0 - bestDistance), Math.round(psHeightCropped / 2.0 - bestDistance), Math.round(2.0 * bestDistance), Math.round(2.0 * bestDistance));
	
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
