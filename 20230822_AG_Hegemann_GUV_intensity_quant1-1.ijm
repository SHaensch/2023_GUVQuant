#@ String (visibility=MESSAGE, value="GUVQuant V1.1 - Sebastian Hänsch" 																					) topMsg
#@ Number (label = "Backgroundwidth (0-inf, 5 default, low=wide, high=narrow)", 														value = "5"			) UserFactor
#@ Number (label = "Width of line measurements (1-inf, 20 default, low=narrow, high=wide)", 											value = "20"		) UserFactorWidth
#@ Number (label = "Maxima detection tolerance (1-inf, 450 default, low=high tolerance, high=low tolerance)", 							value = "450"		) UserTolerance
#@ Number (label = "Maxima detection width (1-inf, default 4, low=wide, high=narrow)", 													value = "4"			) UserMaxWidth
#@ String (visibility=MESSAGE, value="" 																													) empty
#@ String (label = "What is your reference channel", 									choices={"Ch1_blue","Ch2_green","Ch3_red"}, 	style="list"		) UserReference
#@ String (label = "What is your readout channel", 										choices={"Ch1_blue","Ch2_green","Ch3_red"}, 	style="list"		) UserReadout
#@ String (visibility=MESSAGE, value="" 																													) empty
#@ String (label = "Which frames I want to analyze (comma seperated, no spaces, type ´all´ for all frames)",							value = "1,11,21"	) UserFramesAnalysis
#@ String (visibility=MESSAGE, value="" 																													) empty
#@ Boolean(label = "I want to register channels based on signal?"																							) UserChanRegister
#@ Number (label = "Median radius for register blur (0-inf, default 5 for low binding protein, 0 turns it off for clear binding)", 	value = "5"				) UserRegMedian
#@ Number (label = "Variance radius for register blur (0-inf, default 2 for low binding protein, 0 turns it off for clear binding)", 	value = "2"			) UserRegVariance
#@ String (visibility=MESSAGE, value="" 																													) empty
#@ Boolean(label = "I have already created file-associated ROIs that I want to reload?"																		) UserPreviousROILoad
#@ Boolean(label = "I want to automatically save results and store ROI´s for later use?"																	) UserAutosave
#@ String (visibility=MESSAGE, value="" 																													) empty
#@ File (label = "Specify a file to process", style = "file"																								) UserOrigFile
#@ String (visibility=MESSAGE, value="" 																													) empty
#@ Boolean(label = "I want to hide active processes in main loop?"																							) UserBatch
#@ String (visibility=MESSAGE, value="Click OK to Proceed -->"                                                                            					) botMsg   

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* 	Last update: 			21th Aug 2023 V1.1
 *	Main author: 			Sebastian Hänsch, CAI - Center for Advanced Imaging - HHU
 *	Modification author: 	
 *	--- Please consider appropriate acknowledgement and/or coauthorship, if this macro is notably contributing to the success of a publication ---
 *
 *	Description:	Quantification of GUV associated Signal tool
 *					The tool is designed to measure maxima and background of GUV-samples.
 *					The maximum is defined by a reference channel that outlines the GUV strongly and robust.
 *					The found maximum reads the values of the redout channel on this specific maximum position afterwards. 
 *					Background of outer medium and lumen of the GUV is calculated as mean.
 *					ROIs must be set by definition in a first run and are stored for later use / reuse. The are stored within file location.
 *					For channel correction Multistackreg-Plugin is used. The readout channel will be aligned to the reference channel. 
 *						For weak binding protein measurements a median and variance filtering is recommended. For registration of good binding proteins no median or variance is recommended, 
 *						since it could even cause a small mismatch of the channels
 *					
 *	Dependencies: For channel registration: Turbostack-Reg and Multistackreg MUST be installed as plugins (http://bradbusse.net/downloads.html & http://bigwww.epfl.ch/thevenaz/turboreg/)
 *					
 *					Input: 	Any 3 color image would would work. By desing the tool was supposed to have fixed color-shemes associated with it. That is, why the channels are named like that in the macro.
 *							Second Input would be ROI-zips, if they were already created in the first round (see corresponding option). They are named exactly after the file in use.
 *					Output: Creates a results-folder at the file location
 *							1 tiff containing an overview of all ROI´s placed with numbers
 *							2 tiff / stack of the crop region around the first ROI before the correction, ROI to check for registration also over time, if selected by the user
 *							3 tiff / stack of the crop region around the first ROI after the correction, ROI to check for registration also over time, if selected by the user
 *							txt. file with logs for reviewing the parameters used in the first place and bug tracking checkpoints
 *							txt. file of the results from the custom table
 *	
 *	Last modification: 	/SH 1.1 Reused the original version from 2018 and introduced following features:
 *								- added custom output as table
 *								- added choiche of channels and read / reference channels
 *								- added channelregistration with preprocessing options
 *								- added storage and reuse of ROI-sets per file
 *								- added feature to analyze multiple, single or all frames of a timeseries
 *								- added legends for the plots
 *								- added output images in general for monitoring quality of registration, ROI placement, and maxima-detection compare between reference and readout channels
 *								- general cleanup and simplification of the macro
 *								- adding more detailed annotations
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Functions///////////////////////////////////////////////////////////////////////////////////

//-//

//Macro-Script///////////////////////////////////////////////////////////////////////////////

//Clean everything which is still open by chance
while (nImages>0) 
{ 
	selectImage(nImages); 
	close(); 
}
run("Clear Results");
roiManager("reset");
close("Log");

//Create necessary file seperators
FSeperator = File.separator;
print(FSeperator);

//Grab time for timespamping the results
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
month = month+1;
starttime=getTime();

//Open File
run("Bio-Formats", "open=["+UserOrigFile+"] autoscale color_mode=Composite display_rois rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");

//Create Empty Table for later use
PlaceHolder = newArray();
Table.create("Results"); //new V1.1 use name "Results" to set it active
Table.setColumn("ImageName", PlaceHolder);			
Table.setColumn("Frame", PlaceHolder);
Table.setColumn("ROI", PlaceHolder);
Table.setColumn("FrontAVG", PlaceHolder);
Table.setColumn("Maximum", PlaceHolder);
Table.setColumn("BackAVG", PlaceHolder);
Table.setColumn("Fold-Enrich MaxToBack", PlaceHolder);
IJ.renameResults("ResultsTableOfAnalysis"); //new V1.1 use different name to set it inactive

//Create necessary variables for storage and use
SumFront = newArray();
SumBack = newArray();
SumMax = newArray();
ErrorCount = 0;
TableRowCount = 0;
name = getInfo("image.filename");
getPixelSize(unit, pixelWidth, pixelHeight);
TempPath = getDirectory("temp");
nameWOextension = File.getNameWithoutExtension(UserOrigFile);
UserROIPath = File.directory+nameWOextension+".zip";

//Create temporary paths for registration matrix and the resultsfolder
TempPath = getDirectory("temp");
TempFileMatrix = TempPath+"TransformationMatrices.txt";
ResFolder = File.directory+"Results_"+year+"_"+month+"_"+dayOfMonth+"__"+hour+"_"+minute+"_"+second+nameWOextension+FSeperator;

// Create results directory
File.makeDirectory(ResFolder);

//Text string of user input which frames to analyze is converted to array for later and check also if "all" frames are selected
FrameNumbers = newArray();
if (UserFramesAnalysis != "all")
{
	TempUserFramesAnalysis = UserFramesAnalysis;
	while (TempUserFramesAnalysis.length > 0) 
	{
		//use comma as seperator and add individual numbers
		lastIndex = TempUserFramesAnalysis.lastIndexOf(",");
		print(lastIndex);
		TempAdd = TempUserFramesAnalysis.substring(lastIndex+1,TempUserFramesAnalysis.length);
		print(TempAdd);
		FrameNumbers = Array.concat(TempAdd, FrameNumbers);
		if(lastIndex>=0)
		{
			TempUserFramesAnalysis = TempUserFramesAnalysis.substring(0,lastIndex);
		}
		else
		{
			TempUserFramesAnalysis = "";
		}
	}
}
else 
{
	//fill array with all numbers if "all" is selected
	getDimensions(Origwidth, Origheight, Origchannels, Origslices, Origframes);
	for (i = 0; i <= Origframes; i++) 
	{
		FrameNumbers = Array.concat(FrameNumbers, i);
	}
}
Array.print(FrameNumbers);

//If user has choosen channel correction use MultistackReg for possible correction
if (UserChanRegister == true) 
{
	//Create temporary files to do preprocessing for registration
	selectImage(name);
	run("Duplicate...", "title=OriginalBeforeCorrection duplicate");
	selectImage(name);
	rename("Original");
	run("Split Channels");
	selectWindow("C1-Original");
	rename("Ch1_blue");
	selectWindow("C2-Original");
	rename("Ch2_green");
	selectWindow("C3-Original");
	rename("Ch3_red");	
	selectWindow(UserReference);
	run("Duplicate...", "title=RegisterRef duplicate");
	selectWindow(UserReadout);
	print(TempFileMatrix);
	run("Duplicate...", "title=RegisterRead duplicate");
	
	//Do preprocessing on the temporary files according to user parameters
	if(UserRegMedian != 0)
	{
		run("Median...", "radius="+UserRegMedian+" stack");
	}
	if(UserRegVariance != 0)
	{
		run("Variance...", "radius="+UserRegVariance+" stack");
		run("Median...", "radius="+UserRegVariance+" stack");
	}
	run("Tile");
	
	//Do the actual registration on the preprocessed files, save the matrix temporarily and apply the matrix to the original files
	run("MultiStackReg", "stack_1=RegisterRef action_1=[Use as Reference] file_1=[] stack_2=RegisterRead action_2=[Align to First Stack] file_2=["+TempFileMatrix+"] transformation=[Rigid Body] save");
	run("MultiStackReg", "stack_1="+UserReference+" action_1=[Use as Reference] file_1=[] stack_2="+UserReadout+" action_2=[Load Transformation File] file_2=["+TempFileMatrix+"] transformation=[Rigid Body]");	
	
	//Fuse the registered files
	selectWindow("Ch1_blue");
	run("Blue");
	selectWindow("Ch2_green");
	run("Green");
	selectWindow("Ch3_red");
	run("Red");
	run("Merge Channels...", "c1=Ch1_blue c2=Ch2_green c3=Ch3_red create");
	
	//close / delete temporary files 
	rename(name);
	close("RegisterRead");
	close("RegisterRef");
	File.delete(TempFileMatrix);			
}
	
//Give the user the choice to create new ROIs if not loading old ones
if(UserPreviousROILoad == false)
{
	setTool("line");
	selectWindow(name);
	waitForUser(" Before proceeding:\n Make NOW line selections and add them to the ROI Manager (press ctrl+´t´)!\n \n The ROIs should contain the peaks to measure!\n There must be not more than a single maximum, which should be located in the middle!" );
	roiManager("save", UserROIPath);
}
else 
{
	roiManager("open", UserROIPath);
}

//Prepare File to see all ROIs and save / close it
selectWindow(name);
run("Duplicate...", "title=Overview duplicate");
roiManager("Show All");
roiManager("Show All with labels");
run("Flatten");
saveAs("Tiff",ResFolder+"1_ROIplacement.tiff");
selectWindow("1_ROIplacement.tiff");
close();
selectWindow("Overview");
close();

//Prepare files to monitor correction files on first ROI for all frames before and after the correction / save and close it again
if (UserChanRegister == true) 
{
	selectWindow(name);
	roiManager("select", 0);
	getBoundingRect(x, y, width, height);
	makeRectangle(x-UserFactorWidth, y-UserFactorWidth, width+UserFactorWidth*2, height+UserFactorWidth*2);
	run("Duplicate...", "title=AfterCorrectionSmall duplicate");
	saveAs("Tiff", ResFolder+"3_AfterCorrection.tiff");
	close("3_AfterCorrection.tiff");
	
	selectWindow("OriginalBeforeCorrection");
	roiManager("select", 0);
	getBoundingRect(x, y, width, height);
	makeRectangle(x-UserFactorWidth, y-UserFactorWidth, width+UserFactorWidth*2, height+UserFactorWidth*2);
	run("Duplicate...", "title=AfterCorrectionSmall duplicate");
	saveAs("Tiff", ResFolder+"2_BeforeCorrection.tiff");
	close("2_BeforeCorrection.tiff");
	close("OriginalBeforeCorrection");
}

//put in batch mode if the user activated
if (UserBatch)
{
	setBatchMode(true);
}

//Start Main loop
roiManager("show all without labels");
roiManager("show none");
run("Select None");

for (ROIi=0; ROIi<roiManager("count"); ROIi++)
{	
	
	//Create a smaller Crop for active ROI inculde extended boundaries dependend on the linewidth
	selectWindow(name);
	roiManager("select", ROIi);
	getBoundingRect(x, y, width, height);
	makeRectangle(x-UserFactorWidth, y-UserFactorWidth, width+UserFactorWidth*2, height+UserFactorWidth*2);
	run("Duplicate...", "title=Temporary duplicate");
	
	//splitting image 
	run("Split Channels");
	selectWindow("C1-Temporary");
	rename("Ch1_blue");
	selectWindow("C2-Temporary");
	rename("Ch2_green");
	selectWindow("C3-Temporary");
	rename("Ch3_red");	
	
	//starting frame loop for a specific ROI at a time on the crop
	for (Framei = 0; Framei<FrameNumbers.length; Framei++) 
	{
		//generating temporary plot in the reference channel and storing values for later use
		selectImage(UserReference);
		roiManager("select", ROIi);
		Stack.setFrame(FrameNumbers[Framei]);
		Stack.getPosition(channel, slice, frame);
		frameNumber = frame;
		Roi.setStrokeWidth(UserFactorWidth);
		run("Plot Profile");
		rename("RefPlot");
		Plot.getValues(xReference, yReference);
				
		//generating temporary plot in the readout channel and storing values for later use
		selectImage(UserReadout);
		roiManager("select", ROIi);
		Stack.setFrame(FrameNumbers[Framei]);
		Stack.getPosition(channel, slice, frame);
		frameNumber = frame;
		Roi.setStrokeWidth(UserFactorWidth);
		run("Plot Profile");
		rename("ReadPlot");
		Plot.getValues(xReadout, yReadout);
				
		//Finding maxima in reference channel according to user data and defining with saceld widths for front and back end calculations
		length = xReference.length;
		Backgroundlength = parseInt(round(length/UserFactor));
		LeftBorder = parseInt((length/2)-Backgroundlength/UserMaxWidth);
		RightBorder = parseInt((length/2)+Backgroundlength/UserMaxWidth);
		maximumArray= Array.findMaxima(yReference, UserTolerance);
		maximumPos=0;
		noMax = maximumArray.length;
		
		//check if there are maxima at all, if not count as error and return message
		if(noMax == 0)
		{
			print("ERROR!!! Reference ROI"+ROIi+". Has NO Significant Maximum!");
			Array.getStatistics(yReference,min,max,mean,stdDev);
			print("Background: "+mean+", and was not taken for overall statistics!\n  "); 
			ErrorCount = ErrorCount+1;	
		}
		else
		{
			//checking the maxima if they are between the user given limits in the middle. If not count as error.
			for (Maxi=0; Maxi<noMax; Maxi++)
			{
				//Indetify Maxima x and read out the value of the readout channel on that position 
				maximumPos=maximumArray[Maxi];
				maximumY=yReadout[maximumPos];
				maximumX=maximumPos;
				frontEnd = newArray();
				backEnd = newArray();
							
				//If maximum is within the given limits add it to the overall results and calculate Front and back averages		
				if (maximumX >= LeftBorder && maximumX <= RightBorder)
				{
					SumMax = Array.concat(SumMax, maximumY);
					frontEnd = Array.slice(yReadout, 0, Backgroundlength);
					backEnd = Array.slice(yReadout, (length-Backgroundlength), length);
					Array.getStatistics(frontEnd,min,max,mean,stdDev);
					SumFront = Array.concat(SumFront, mean);
					print(mean);
					print(maximumY);
					Array.getStatistics(backEnd,min,max,mean,stdDev);
					print(mean);
					SumBack = Array.concat(SumBack, mean);
					Maxi=noMax;
					MaxFound = 1;
				}
				else
				{
					MaxFound = 0;	
				}
			}
			
			//return error, if there might be a max, but not in the limited center
			if (MaxFound == 0)
			{
				print("ERROR!!! No Maximum of Reference measurement "+ROIi+" within the center limits!");
				Array.getStatistics(yReadout,min,max,mean,stdDev);
				print("Background: "+mean+", and was not taken for overall statistics!\n  ");
				ErrorCount = ErrorCount+1;
			}
		}
					
		//generating the actual plot for displaying starting with the plot of the output channel
		nametif = replace(name, ".nd2", "");
		selectWindow("ReadPlot");
		close();
		selectWindow("RefPlot");
		close();
		Plot.create("Measurement "+UserReadout+" "+ROIi, unit, "Grey values", xReadout, yReadout);
		Plot.setColor("black");	
		
		// Is there a maximum or not?
		if(noMax == 0)
		{
			//leave out the plots if no maxima were found and calculating overall background instead
			Plot.setColor("red");
			AbsBackLength = Backgroundlength*pixelWidth;
			yValues = newArray(mean, mean); 
			xValues = newArray(0, xReadout.length);
			Plot.add("line", xValues, yValues);
		}
		else
		{
			//draw line for the reference line profile in grey
			Plot.setColor("gray");
			yValues = yReference; 
			xValues = xReference;
			Plot.add("line", xValues, yValues);
			
			//draw front end background horizontal line in red
			Plot.setColor("red");
			AbsBackLength = Backgroundlength*pixelWidth;
			Array.getStatistics(frontEnd,min,max,mean,stdDev);
			yValues = newArray(mean, mean); 
			xValues = newArray(0, AbsBackLength);
			Plot.add("line", xValues, yValues);
			
			//draw front end vertical line
			xValues = newArray(AbsBackLength, AbsBackLength);
			yValues = newArray(0, mean); 
			Plot.add("line", xValues, yValues);
			
			//draw back end horizontal line
			AbsBackLength = (length*pixelWidth)-AbsBackLength;
			Array.getStatistics(backEnd,min,max,mean,stdDev);
			yValues = newArray(mean, mean); 
			xValues = newArray(AbsBackLength, length);
			Plot.add("line", xValues, yValues);
			
			//draw back end vertical line
			xValues = newArray(AbsBackLength, AbsBackLength);
			yValues = newArray(0, mean); 
			Plot.add("line", xValues, yValues);
			
			//draw maximum with readout height if there is one maximum, otherwise draw baseline
			AbsBackLength = maximumX*pixelWidth;
			height = yReadout[maximumX];
			yValues = newArray(0, height); 
			xValues = newArray(AbsBackLength, AbsBackLength);
			Plot.setColor("blue");
			if (MaxFound != 0)
			{	
				Plot.add("line", xValues, yValues);
			}
			else
			{
				
			}
			
			//draw left border of maximum detection
			AbsBackLength = LeftBorder*pixelWidth;
			yValues = newArray(0, 255); 
			xValues = newArray(AbsBackLength, AbsBackLength);
			Plot.setColor("blue");
			Plot.add("line", xValues, yValues);
			
			//draw right border of maximum detection
			AbsBackLength = RightBorder*pixelWidth;
			xValues = newArray(AbsBackLength, AbsBackLength);
			Plot.setColor("blue");
			Plot.add("line", xValues, yValues);
				
		}
		
		//Increase resolution a bit to see plots better
		Plot.setColor("black");
		Plot.show();
		Plot.makeHighResolution("Temporary Plot",1.5);
		
		//insert plot information strings and flatten
		setFont("Calibri", 12);
		drawString("Frame: "+Framei, 0, 12, "black");
		drawString("ROI: "+ROIi, 0, 24, "black");
		drawString("Black: "+UserReadout, 0, 36, "black");
		drawString("Grey: "+UserReference, 0, 48, "black");
		run("Flatten", "");
	
		//concatenate plots to a single multi tiff file output
		if (TableRowCount==0)
		{
			rename("Plots "+UserReadout+" "+ROIi);
			print("first image");
			selectImage("Temporary Plot");
			close();
		}
		else 
		{
			rename("Temporary Plot2");
			run("Copy");
			selectImage("Plots "+UserReadout+" 0");
			run("Add Slice");
			run("Paste");
			selectWindow("Temporary Plot");
			close();
			selectImage("Temporary Plot2");
			close();
		}
		
		//calculate EnrichmentScore
		MaxTOBackEnrich = SumMax[TableRowCount]/SumBack[TableRowCount];
		
		//fill in results for the custom table and sort it by ROI to track changes over time better
		selectWindow("ResultsTableOfAnalysis");
		IJ.renameResults("Results"); //activate table
		setResult("ImageName", TableRowCount, name);
		setResult("Frame", TableRowCount, frameNumber);
		setResult("ROI", TableRowCount, ROIi);
		setResult("FrontAVG", TableRowCount, SumFront[TableRowCount]);
		setResult("Maximum", TableRowCount, SumMax[TableRowCount]);
		setResult("BackAVG", TableRowCount, SumBack[TableRowCount]);
		setResult("Fold-Enrich MaxToBack", TableRowCount, MaxTOBackEnrich);
		Table.sort("ROI");
		updateResults();
		IJ.renameResults("ResultsTableOfAnalysis"); //deactivate table
		
		//Checkpoint for Logs
		print("\n\n "+UserReadout+" Result for measurement frame:"+Framei+" ROI: "+ROIi+"\n MaxGray: "+SumMax[TableRowCount]+"\n BAckground at Front: "+SumFront[TableRowCount]+"\n BAckground at Back: "+SumBack[TableRowCount]);
		
		//Increase Tablecount to have an incrementing variable for result placement for both frames together with ROIs
		TableRowCount++;
		
		//Close leftover of the plot
		close("Measurement "+UserReadout+" "+ROIi);
	}
	
	//close all crops for the next round and next ROI
	selectWindow("Ch1_blue");
	close();
	selectWindow("Ch2_green");
	close();
	selectWindow("Ch3_red");
	close();
		
}		

//returning final calculation for the different parameter and put to log...does not make sense if frames are selected, so deactivated 
/*
AbsBackLength = Backgroundlength*pixelWidth;
print("\n ____________________________________________________________________ \n Overall Statistics"+UserReadout+": ");
Array.getStatistics(SumMax,min,max,mean,stdDev);
print("Maxima:\n Mean: "+mean+"   Min: "+min+"   Max: "+max+"   stdDev: "+stdDev);
Array.getStatistics(SumFront,min,max,mean,stdDev);
print("Front Background :\n Mean: "+mean+"   Min: "+min+"   Max: "+max+"   stdDev: "+stdDev);	
Array.getStatistics(SumBack,min,max,mean,stdDev);
print("Back Background :\n Mean: "+mean+"   Min: "+min+"   Max: "+max+"   stdDev: "+stdDev);
print("Calculated by "+SumMax.length+" measurements. "+ErrorCount+" measurements did not show a maximum in the middle and counted as Error!");	
*/

//Document user input parameters
print("__________________________________________________");
print("Used parameters:");
print("Backgroundwidth: "+UserFactor);
print("Width of line measurements: "+UserFactorWidth);
print("Maxima detection tolerance?: "+UserTolerance);
print("Maxima detection width?: "+UserMaxWidth);
print("What is your reference channel: "+UserReference);
print("What is your readout channel: "+UserReadout);
print("Which frames I want to analyze: "+UserFramesAnalysis);
print("I want to register channels based on signal?: "+UserChanRegister);
print("I have already created file-associated ROIs that I want to reload?: "+UserPreviousROILoad);
print("I want to automatically save results and store ROI´s for later use?: "+UserAutosave);
print("Specify a file to process: "+UserOrigFile);
print("__________________________________________________");
print("Results saved to:");
print(ResFolder);

//Calculate final time, in order to track performance and possibilities to speed up
endtime=getTime();
timeDifference=endtime-starttime;
seconds = timeDifference/1000;
minutes = floor(seconds/60);
seconds = seconds - minutes*60;
hours = floor(minutes/60);
minutes = minutes - hours*60;
print("\nIt took "+hours+" hours, "+minutes+" minutes and "+seconds+" seconds for calculation!");

//Save all files
selectWindow("Log");
save(ResFolder+"ParamenterLog.txt");
selectWindow("ResultsTableOfAnalysis");
saveAs("Text", ResFolder+"ResultsFor"+nameWOextension);
selectImage("Plots "+UserReadout+" 0");
saveAs("Tiff", ResFolder+"4_"+nameWOextension+"_Plots.tiff");
print("\\Clear");
close("Log");

//exit batchmode if activated by user
if (UserBatch)
{
	setBatchMode(false);
}

//Clean everything up 
close("*");

//Finish
beep();
showMessage("Analysis Done!");
	