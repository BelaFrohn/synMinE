// Min_wave_analysis_in_E.coli
// ImageJ macro for detection of E.coli cells with fluorescently labeled Min oscillation, make kymographs, and analyze the oscillation periods
// Validated only in Fiji v1.53f or later
// Input: 2D time-lapse image files with fluorescently labeled Min oscillation inside E.coli cells
// Resultant data sets are saved as new Tiff image files (for kymographs), a CSV file (for oscillation period analysis), and/or a Zip file (for ROIs)
// For detailed information: Kohyama et al, Nature Communications 2024
//
// Shunshi Kohyama (Max Planck Institute of Biochemistry), 2024

// Presetting
methods = getList("threshold.methods");
getVoxelSize(width, height, depth, unit);
interval =  Stack.getFrameInterval();
Dialog.create("Min wave analysis");
Dialog.addString("Title", getTitle());
if (unit != "µm") {
	Dialog.addMessage("                          Please check!!", 12, "red");
}
Dialog.addNumber("Set scale (µm in a pixel)", width);
if (interval == 0) {
	Dialog.addMessage("                          Please check!!", 12, "red");
}
Dialog.addNumber("Set interval (sec)", interval);
Dialog.addNumber("Set channel (color)", 1);
Dialog.addChoice("Threshold method", methods, "Huang");
labels = newArray("Wave analysis", "Save results (csv)", "Save kymographs", "Save ROIs");
defaults = newArray(true, true, true, true);
answer = newArray(4);
Dialog.addCheckboxGroup(2, 2, labels, defaults);
Dialog.show();
title = Dialog.getString();
scale = Dialog.getNumber();
interval = Dialog.getNumber();
channel = Dialog.getNumber();
type = Dialog.getChoice();
for (i=0; i<4; i++) {
	answer[i] = Dialog.getCheckbox();
}
	
if (answer[0] == 1) {
	smoothingArray = newArray(5, 10, 20, 40);
	Dialog.create("Oscillation analysis");
	Dialog.addRadioButtonGroup("Smoothing parameter", smoothingArray, 2, 2, "10.0");
	Dialog.addNumber("Max. cell width (µm)", 2);
	Dialog.addNumber("Max. cell length (µm)", 5);
	Dialog.addNumber("Max. period (s)", 60);
	Dialog.show();
	smoothing = Dialog.getRadioButton();
	maxWidth = Dialog.getNumber();
	maxlength = Dialog.getNumber();
	maxPeriod = Dialog.getNumber();
}

if (answer[1] == 1 || answer[2] == 1 || answer[3] == 1) {
	dir = getDirectory("Choose a Directory");
}

if (isOpen("ROI Manager")) {
selectWindow("ROI Manager");
run("Close");
}
run("Set Scale...", "distance=1 known=scale pixel=1 unit=µm");

setBatchMode(true);

// Cell detection
run("Set Measurements...", "  redirect=None decimal=3");
run("Duplicate...", "duplicate channels=channel");
run("Z Project...", "projection=[Average Intensity]");
setAutoThreshold(type+" dark");
run("Analyze Particles...", "size=0.1-Infinity exclude clear add");
close(); 
n = nSlices;
N = roiManager("count");
run("Clear Results");

// Oscillation analysis
for (i=0; i<N; i++) {
	showProgress(-i/N);
	roiManager("select", i);
	angle = getValue("Angle");
	run("Duplicate...", "duplicate");
	Roi.remove;
	if (angle > 90) {
		rotation = angle - 180;
	} else {
		rotation = angle;
	}	
	run("Rotate... ", "angle=rotation grid=1 interpolation=Bilinear enlarge stack");
	run("Z Project...", "projection=[Average Intensity]");	
	setAutoThreshold(type+" dark");
	run("Analyze Particles...", "add");
	close();
	roiManager("select", N);
	majorAxis = getValue("Major");
	minorAxis = getValue("Minor");
	setResult("ID", i, title+"_"+i+1);
	setResult("Major axis (µm)", i, majorAxis);
	setResult("Minor axis (µm)", i, minorAxis);
	
	// Kymograph
	if (answer[0] == 1 || answer[2] == 1) { 
		run("Duplicate...", "duplicate");
		w = getWidth();
		run("Size...", "width=w height=1 depth=n average interpolation=Bilinear");
		run("Make Montage...", "columns=1 rows=n scale=1");
		if (answer[2] == 1) {
		resetMinAndMax();
		saveAs("Tiff", dir+title+"_"+i+1+".tif");
		}
	}

	// Oscillation period detection
	if (answer[0] == 1) {
		width = w * smoothing;
		run("Size...", "width=width height=n depth=1 average interpolation=Bilinear");
		int = 0;
		pos = 0;
		makeRectangle(0, 0, 5, n);
		for (j = 0; j <= w; j++) {
			Roi.move(j, 0);
			temp = getValue("Mean");
			if (temp > int) {
				int = temp;
				pos = j;
			}
		}
		makeRectangle(pos, 0, 5, n);
		run("Duplicate...", " ");
		height = n * smoothing;
		run("Size...", "width=1 height=height depth=1 average interpolation=Bilinear");
		makeLine(0, 0, 0, height);
		profile = getProfile();
		close();
		
		timepoints = Array.getSequence(lengthOf(profile));
		for (j=0; j<lengthOf(profile); j++) {
			timepoints[j] = j * interval / smoothing;
		}
		Fit.doFit("y = a * sin(b * x + c) + d", timepoints, profile);
		R = Fit.rSquared;
		periodRaw = Fit.p(1);
		repeat = Math.floor(periodRaw/PI);
		periodAmp = periodRaw - PI * repeat;
		if (periodAmp > PI/2) {
			periodAmp = PI - periodAmp;
		}
		period = 2 * PI / periodAmp;
		setResult("Period (sec)", i, period);
		setResult("Fitting (R^2)", i, R);
		if (majorAxis > maxlength || minorAxis > maxWidth || period > maxPeriod) {
			setResult("Detection", i, "Check!!!");
		} else {
			setResult("Detection", i, "");
		}
	}
	
	if (answer[0] == 1 || answer[2] == 1) {
		close();
		close();
	}
	close();
	
	while (roiManager("count") > N) {
		roiManager("select", N);
		roiManager("delete");
	}

}

if (answer[3] == 1) {
	roiManager("save", dir+title+"_ROI.zip")
}
if (answer[1] == 1) {
	saveAs("Results", dir+title+".csv");
}

close();
setBatchMode(false);