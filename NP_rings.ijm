function findRoisWithName(roiName) { 
	nR = roiManager("Count"); 
	roiIdx = newArray(nR); 
	k=0; 
	clippedIdx = newArray(0); 
	
	for (i=0; i<nR; i++) { 
		roiManager("Select", i); 
		rName = Roi.getName(); 
		if (startsWith(rName, roiName)) { 
			roiIdx[k] = i; 
			k++; 
		} 
	} 
	if (k>0) { 
		clippedIdx = Array.trim(roiIdx,k); 
	} 
	 
	return clippedIdx; 
} 

SourceDir = getDirectory('~/');
list = getFileList(SourceDir+"img"+File.separator);
File.makeDirectory(SourceDir+"mask_ab");
File.makeDirectory(SourceDir+"mask_NP");
File.makeDirectory(SourceDir+"mask_nonNP");
File.makeDirectory(SourceDir+"ab_summary");

for(j = 0; j < list.length; j++) {
	open(SourceDir+"img"+File.separator+list[j]);
	title = getTitle();
	AddChannel(title);
}

function AddChannel(title) { 

//title = getTitle();
run("Set Scale...", "distance=1 known=1 unit=micrometer global");
run("Set Measurements...", "area mean min redirect=None decimal=3");

//segment 4G8 channel 
setSlice(3);
run("Duplicate...", "  channels=3");
rename("4G8");
run("Enhance Contrast", "saturated=0.35");
run("Despeckle", "stack");
setOption("ScaleConversions", true);
run("Gaussian Blur...", "sigma=5");
setOption("BlackBackground", false);
setAutoThreshold("Otsu dark");
//run("Analyze Particles...", "size=50-Infinity pixel show=Masks display summarize add slice");
run("Analyze Particles...", "size=100-Infinity pixel show=Masks display add");
	//selectWindow("Summary of 4G8");
	//Table.set("Slice", 0, "4G8"); 
selectWindow("Mask of 4G8");
run("Create Selection");
if ( selectionType() != -1) {
Amyloid_Status = "present";
run("ROI Manager...");
roiManager("Add");
n=RoiManager.size;
roiManager("Select", n-1);
roiManager("Rename", "4G8");
//close("Mask of "+title);
//roiManager("Select", 0);
//run("Enlarge...", "enlarge=radius pixel");
//roiManager("Add");
//roiManager("Select", 1);
//roiManager("Rename", region_of_interest+radius+"microns");
}
else {
Amyloid_Status = "empty";
}
selectWindow("Mask of 4G8");
run("Invert LUTs");
saveAs("tiff", SourceDir+"mask_ab"+File.separator+ title);

//segment AT8 channel 
selectWindow(title);
run("Select None");
setSlice(23);
run("Duplicate...", "  channels=23");
rename("AT8");
run("Enhance Contrast", "saturated=0.35");
run("Despeckle", "stack");
setOption("ScaleConversions", true);
setOption("BlackBackground", false);
setAutoThreshold("Otsu dark");
run("Analyze Particles...", "pixel show=Masks summarize slice");
	//selectWindow("Summary of AT8");
	//Table.set("Slice", 1, "AT8"); 
selectWindow("Mask of AT8");
run("Create Selection");
if ( selectionType() != -1) {
Tau_Status = "present";
run("ROI Manager...");
roiManager("Add");
n=RoiManager.size;
roiManager("Select", n-1);
roiManager("Rename", "AT8");
}
else {
Tau_Status = "empty";
}
//close("Mask of "+title);

if (Amyloid_Status != "empty") {
	if (Tau_Status != "empty") {
		//create ROIs of 4G8+ plaque +/- AT8 
		n=RoiManager.size;
		AT8_array=newArray();
		for (i = 0; i < n-2; i++) {
		roiManager("Select", newArray(i,n-1));
		roiManager("AND");
			if (selectionType() != -1) {
				roiManager("Select", i);
				roiManager("Rename", "4G8_AT8"+"_"+i);
				AT8_array = Array.concat(AT8_array,"1");
			}
			else {
				roiManager("Select", i);
				roiManager("Rename", "4G8_NoAT8"+"_"+i);
				AT8_array = Array.concat(AT8_array,"0");
			}
			}
		NP_list = findRoisWithName("4G8_AT8");
		P_list = findRoisWithName("4G8_NoAT8");
		Table.setColumn("AT8", AT8_array,"Results");
		Table.update;
		
	}
}

//create a new image using 4G8 + AT8 if avilable, if not create an empty channel
if (Amyloid_Status != "empty") {
	if (Tau_Status != "empty") {
		roiManager("Select", NP_list);
		roiManager("Combine");
		if ( selectionType() != -1) {
			roiManager("Add");
		}
		n=RoiManager.size;
		roiManager("Select", n-1);
		roiManager("Rename", "4G8_AT8");
		newImage(title, "16-bit black", 500, 500, 1);
		setForegroundColor(255, 255, 255);
		roiManager("Select", n-1);
		run("Fill", "slice");
		saveAs("tiff", SourceDir+"mask_NP"+File.separator+ title);
		//waitForUser("check if image is saved in mask_NP");
	}
	else {
	newImage(title, "16-bit black", 500, 500, 1);
	saveAs("tiff", SourceDir+"mask_NP"+File.separator+ title);
	}
}
else {
	newImage(title, "16-bit black", 500, 500, 1);
	saveAs("tiff", SourceDir+"mask_NP"+File.separator+ title);
}

//create a new image using 4G8 without AT8 if avilable, if not create an empty channel
if (Amyloid_Status != "empty") {
	if (Tau_Status != "empty") {
		roiManager("Select", P_list);
		roiManager("Combine");
		if ( selectionType() != -1) {
			roiManager("Add");
		}
		n=RoiManager.size;
		roiManager("Select", n-1);
		roiManager("Rename", "4G8_NoAT8");
		newImage(title, "16-bit black", 500, 500, 1);
		setForegroundColor(255, 255, 255);
		roiManager("Select", n-1);
		run("Fill", "slice");
		saveAs("tiff", SourceDir+"mask_nonNP"+File.separator+ title);
		//waitForUser("check if image is saved in mask_nonNP");
		
		run("Close All");
		roiManager("reset");
		selectWindow("Results");
		saveAs("Results", SourceDir+"ab_summary"+File.separator+title+".csv");
		close("Results");
		close("Summary");
	}
	else {
		newImage(title, "16-bit black", 500, 500, 1);
		n=RoiManager.size;
		roiManager("Select", n-2); 
		setForegroundColor(255, 255, 255);
		run("Fill", "slice");
		saveAs("tiff", SourceDir+"mask_nonNP"+File.separator+ title);
		run("Close All");
		roiManager("reset");
		selectWindow("Results");
		saveAs("Results", SourceDir+"ab_summary"+File.separator+title+".csv");
		close("Results");
		close("Summary");
	}
}
else {
newImage(title, "16-bit black", 500, 500, 1);
saveAs("tiff", SourceDir+"mask_nonNP"+File.separator+ title);
run("Close All");
roiManager("reset");
//selectWindow("Results");
//saveAs("Results", SourceDir+"ab_summary"+File.separator+title+".csv");
//close("Results");
close("Summary");
}

} //end of function


/* future update
run("Enlarge...", "enlarge=10");
run("Enlarge...", "enlarge=10");
roiManager("Select", newArray(29,31));
roiManager("XOR");
roiManager("Add");
roiManager("Select", 31);
roiManager("Select", 32);
roiManager("Rename", "4G8_AT8_10-20");
roiManager("Select", 31);
run("Enlarge...", "enlarge=10");
roiManager("Add");
roiManager("Select", 33);
roiManager("Rename", "4G8_AT8_30");
*/
