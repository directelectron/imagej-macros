print("\\Clear");
print("DE_DQECurve.ijm");
print("Version : 1.00");
print("Date : 2023.04.03");
print("Author : Benjamin Bammes, Direct Electron LP (bbammes@directelectron.com)");
print("License : GNU General Public License v2.0");
print("");

if (!isOpen("MTF")) {
	print("Error");
	exit("You must have the CSV results of ASI_MTF open in a window titled 'MTF'.");
}
selectWindow("MTF");
mtfFreq = Table.getColumn("Spatial_frequency_(Nyquist)");
mtfData = Table.getColumn("fit");

nnpsExists = false;
if (isOpen("NNPS")) {
	selectWindow("NNPS");
	nnpsFreq = Table.getColumn("X_Mean");
	nnpsData = Table.getColumn("fit");
	nnpsExists = true;
}

dqeZero = getNumber("DQE(0)", 0);
if (dqeZero <= 0) {
	print("Error");
	exit("The value of DQE(0) must be positive.");
}

superRes = getBoolean("Type of data used for MTF calculation", "Super-Res", "Normal Res");

if (superRes) {
	for (i = 0; i < mtfFreq.length - 1; i++) {
		mtfFreq[i] *= 2.0;
	}
}

mtfData[0] = 1.0;

s = 0.0;
sInc = 0.01;
sMax = 1.0;
if (superRes) {
	sMax = 2.0;
}

freq = newArray(0);
mtf = newArray(0);
dqe = newArray(0);
while (s <= sMax) {
	bestIndex = mtfData.length - 1;
	for (i = 1; i < mtfData.length; i++) {
		if (s <= mtfFreq[i]) {
			bestIndex = i;
			break;
		}
	}
	mtfValue = mtfData[bestIndex] + (s - mtfFreq[bestIndex]) * (mtfData[bestIndex - 1] - mtfData[bestIndex]) / (mtfFreq[bestIndex - 1] - mtfFreq[bestIndex]);
	nnpsValue = 1.0;
	if (nnpsExists) {
		bestIndex = nnpsData.length - 1;
		for (i = 1; i < nnpsData.length; i++) {
			if (s <= nnpsFreq[i]) {
				bestIndex = i;
				break;
			}
		}
		nnpsValue = nnpsData[bestIndex] + (s - nnpsFreq[bestIndex]) * (nnpsData[bestIndex - 1] - nnpsData[bestIndex]) / (nnpsFreq[bestIndex - 1] - nnpsFreq[bestIndex]);
	}
	freq = Array.concat(freq, s);
	mtf = Array.concat(mtf, mtfValue);
	dqe = Array.concat(dqe, dqeZero * mtfValue * mtfValue / nnpsValue);
	s += sInc;
}

Plot.create("DQE Results", "1 / Physical Nyquist", "");
Plot.setFrameSize(800, 600);
Plot.setLimits(0, sMax, 0, 1.1 );
Plot.setColor("#bc5090");
Plot.add("line", freq, mtf, "MTF");
Plot.setColor("#003f5c");
Plot.add("line", freq, dqe, "DQE");
Plot.addLegend("MTF\nDQE");
Plot.show();

print("");
print("Done");
