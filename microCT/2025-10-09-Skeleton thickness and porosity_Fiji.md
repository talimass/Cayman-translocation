#Skeleton thickness and porosity_Fiji

##Orientation of the data set

Open data and reslice: image->stack->reslice 

![image](/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Talimass_Lab_Notebook_Mass_Lab/images/Fiji_Reslice.png)

Reslice from the top or from right to make the stack horizontal

![image](/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Talimass_Lab_Notebook_Mass_Lab/images/reslice_stack.png)

Save Substack by using ROI of intrast. ROI can be saved to ROImanager and reopened from there: Analyze->Tools->ROImanager

![image](/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Talimass_Lab_Notebook_Mass_Lab/images/ROImanager.png)

##Data segmentation using trainable weka segmentation

In Fiji, open Plugin -> Segmentation -> trainable weka segmentation

![image](/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Talimass_Lab_Notebook_Mass_Lab/images/Waka.png)

After open rename the labels of the class in setting.


![image](/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Talimass_Lab_Notebook_Mass_Lab/images/weka_open.png)

Using the freehand tool, mark the area at each section and add to the relevant class. Repeat in about 20 slices.

![image](/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Talimass_Lab_Notebook_Mass_Lab/images/weka_seg.png)

Once the segmentation is complete, use the train classifier option.
It is possible to save data, save classifiers, and reload them in the next set.

![image](/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Talimass_Lab_Notebook_Mass_Lab/images/weka_after.png)

Make sure the segmentation correct before creating the results.

After creating the result, adjust the threshold to BW using the Otso. Make sure the void value is 0 and the skeleton is 255, if not, invert the image: edit -> invert.

##Thickness calculation using Bone J

Reslice the binary stack from the bottom to get stack in the same orientation, and that slice one will be the top of the coral. 

Open Fiji and make sure you have the BoneJ plugin; if not, add it. 

###Fiji
Open binary data set the scale at: analyze->set scale in micron

Open: Plugin->bonej->slicegeometry and set the parameters: 
![image](/Users/talimass/Documents/Documents - MacBook Pro/GitHub/Talimass_Lab_Notebook_Mass_Lab/images/bonej setting.png)

Plot the Mean Thick 3D (¬µm) by slice from the csv file

##Thickness and porosity using Python

On Coral1: open command prompt (CMD) 
Run: C:\ProgramData\miniconda3\Scripts\activate.bat
Create a working conda environment: as instructed 
[here] (https://youtu.be/hDGSZMLS5F4?t=342)

After creating a working environment, install Python on it 	as instructed 
[here] (https://youtu.be/hDGSZMLS5F4?t=342)

Activate the environment by: conda activate XXX

Copy the path of the Python code and run it






