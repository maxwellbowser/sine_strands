# sine_strands 
### Description:
Python/MATLAB Analysis pipeline for Sinusoidal Perturbation.  
Extracts elastic moduli, viscus moduli, and individual components of complex moduli (A, k, B, C, and rates) from sinusoidal perturbation data.

### Installation:
**Python Files**  
For compiled version, download the zip file containing all files, then then create a shortcut on your desktop which launches the sine_strands application (.exe).
It is important that the application file DOES NOT leave the folder, because otherwise it won't know where to find the software files. I reccomend having the compiled folder somewhere accessible (since all processed data will be saved in this file) but safe, so you don't accidentally add any files into the folder (which can break the software). 

For running Python files, ensure that all used libraries are installed and then run the script using VScode or other Python interface.  

**MATLAB Files**  
This MATLAB software is only tested on version R2023b, so ensure that MATLAB R2023b is installed.
MATLAB files are included in a folder, which has all of the needed MATLAB files included, as well as an output folder. Download this folder and place it somewhere accesible. 

## File Naming Convention:  
To begin, sine_strands uses a very strict naming convention to extract file data properly.  
This naming convention is described here, where options for words are included in brackets:  

**SinusoidalPerturbation[WT/KO][Mouse #][Treated/Untreated][Relax/Active %][-SL:optional]**  
          Only include Sarcomere Length (SL) if NOT 225

Examples: 
--
WT Mouse ID # 1, Treated, Relaxed, @ 225 SL  
**SinusoidalPerturbationWT01TreatedRelax.dat**

KO Mouse ID # 4, Untreated, Active @ 86%, @ 200 SL  
**SinusoidalPerturbationKO04Untreated86-200.dat**

WT Mouse ID # 988, Untreated, Active @ 72%, @ 225 SL  
**SinusoidalPerturbationWT988Untreated72.dat**

Please note that "SinusoidalPerturbation" is NOT required to be present, and can be swapped with any other text (so long as it does not include any of the other described file information)  

  
For example:  
WT Mouse ID # 488, Untreated, Active @ 76%, @ 200 SL  

could have the following filenames:  
1) SinusoidalPerturbationWT488Untreated76-200 
2) SinePertWT488Untreated76-200   
3) WT488Untreated76-200  
4) **etc.**

## Sine_Strands Usage:  
Using the compiled version of sine_strands, open the application "sine_strands" or double-click a shortcut linked to the "sine_strands" application file, located with the compiled files. (Found in "Releases" in zipped file)

This will open the GUI:  
<img src="https://github.com/user-attachments/assets/e8a43ce0-8f8f-4cdf-87a7-07fd0a8db99c" alt="SineStrands GUI" style="width:25%; height:auto;">

#### Simple (Uncalibrated):  
If you want to simply analyze a sinusoidal perturbation file, enter your desired output folder name (e.g. "Example_Sample-2025_01_01") in the "Folder Name:" box then click **"Open Files"**.
A file selection screen will open up, where you can select any number of .DAT files for analysis, and click "Open". A few moments later, Sine_strands should close and a folder with your analyzed files will open:  
  
<img src="https://github.com/user-attachments/assets/d5c20fe4-ee67-4b24-b6b1-15ffe1fc96e0" alt="SineStrands GUI" style="width:50%; height:auto;">  
  

This folder contains a summary excel doc of all analyzed files (if multiple files were analyzed) elastic/viscous moduli:  
  
<img src="https://github.com/user-attachments/assets/de723bc5-3a31-4a3f-a3bf-b1962c93895f" style="width:60%; height:auto;">  

    

This file is used for the downstream analysis using MATLAB. To move onto the MATLAB anaysis portion, please continue to "MATLAB Analysis Usage"

The Sine_strands output folder also contains a folder with all individual files' in-depth analysis data as .csv files:  
  
<img src="https://github.com/user-attachments/assets/1416d01a-1d10-456e-8217-6ef3092e26e6" style="width:60%; height:auto;">  

Please note that the "Uncalibrated" Vm/Em columns are identical to the other Em/Vm columns in the above example. This is because a calibration file was not used for this analysis!

#### Calibrated:  
If you have a calibration file you'd like to use (which will subtract off any moduli, removing delay), you will need to first run the calibration file.
Therefore, follow the above instructions for "simple" analysis using the calibration DAT recording and save the in-depth .CSV file (2nd Spreadsheet shown above). The sine_strands program will ONLY accept this .CSV file for calibration, to prevent users from accidentally inputting the summary file (which is .xlsx) for calibration. Sine_strands REQUIRES the in-depth data to properly calibrate sinusoidal input data.  

**Note:** Unless you added some sizes, you should get a little error pop-up which tell you you're missing size measurements for your cell. This just means sine_strands couldn't find dimensions for your cell, which is good because calibration files should not be run on cells.  

After you have analyzed the calibration file, you will re-open Sinestrands and click "Select" underneath "Calibration File" in the GUI. This will open a file selection window, where you can navigate to the calibration file and select it, then hit "Open". (NOTE: You may get an error about "calibration file" not being found, please read the "IMPORTANT" section below)
The file navigation window should now close, and you'll be presented with an updated GUI, which lists the calibration name under "Calibration File:"  
**GUI Before selecting calibration file:**  
  
<img src="https://github.com/user-attachments/assets/e8a43ce0-8f8f-4cdf-87a7-07fd0a8db99c" alt="SineStrands GUI" style="width:25%; height:auto;">

**GUI After selecting calibration file:**  
  
<img src="https://github.com/user-attachments/assets/eb13cdf4-be19-47b7-b819-330708b6f379" alt="SineStrands GUI" style="width:25%; height:auto;">  

**IMPORTANT - FOR THE FIRST TIME ONLY:**  
This selected calibration file is now selected, but the program must be reset before the file can be used. Therefore, close out sine_strands and relaunch. Now when you relaunch sine_strands, you should have the calibration file displayed. Sine_strands will remember this calibration file, so long as it does not change locations.
  
Now you can follow the rest of the example above, where you'll put in a folder name and hit "Select Files", then select the files you'd like to analyze. **Please note that the calibrated file must have IDENTICAL frequency sweeps to the analyzed files.**  



Once you have your calibrated summary file, you can continue onto the MATLAB analysis!  
  

## MATLAB Analysis Usage:  

To preface this section, the MATLAB analysis requires more hands-on work with code files. BEWARE!! 👻  

1) Start by copying/moving your output summary file from Sine_strands "Data_Summary-ExampleName" into your folder containing the MATLAB files.
2) Open "analyze_Nyquist_plots" in MATLAB, which should give a view like this:
<br>

![image](https://github.com/user-attachments/assets/49a5319a-3f07-437a-b9cb-c5cb04935587)  

<br>
  
3) Change the input_dir text on line #7 between the single quote marks (') to the path of your MATLAB folder where this "analyze_Nyquist_plots" file is located.  (You only have to do this once, unless you move this whole folder!!)

4) Now, copy the name of your output file into the single quotes on line #10, changing filename to the name of your summary file from SineStrands (e.g. filename = 'Data_Summary-ExampleName.xlsx')
4a) BE SURE TO ADD ".xlsx" to the end of the filename!

5) Change "heart_sample" on line #18 to the mouse ID that you want to extract data from. The output from sinestrands always adds "Mouse" before the number, i.e. "Mouse04" or "Mouse988" would be the heart sample.

6) Change to the "EDITOR" tab at the top, if you're not already there, then click "Run". This should output a few graphs and then automatically stop.
<br>
   
![image](https://github.com/user-attachments/assets/dcde2bba-0add-4354-881d-c8f527e2c11a)  


<br>


8) For this sample, the output data will be created in the base folder, and the graphs will be saved to the "output" folder that is inside this folder.

**Output excel file is created in base folder, graphs are placed in "output" folder:**  
<br>


<img src="https://github.com/user-attachments/assets/313e1224-e0bf-4bf3-b53c-dac0d07d4187" style="width:50%; height:auto;">   

**Output data from MATLAB analysis:**  

![image](https://github.com/user-attachments/assets/901a64bb-12b2-41dd-ad28-f1a07e69e827)  


10)   To run the rest of the samples in the sinestrands file, move the output file to a new place (or it might get overwritten), then change "heart_sample" to the next mouse ID # and hit run again.


**Important:**  
For the next file you analyze, you only need to change the name of the summary file on line #10, provided that you've copied the summary file into the folder with the MATLAB files.


And that's it! From here, you can use the different curve fit metrics of the nyquest plot to describe how two tissues behave differently to mechancial forces.





