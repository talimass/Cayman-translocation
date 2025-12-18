# -*- coding: utf-8 -*-

"""
#####################################################################
# Developed in 2025 by Isabela Vitienes, in the Zaslansky lab 
# of the Department of Operative, Preventive and Pediatric Dentistry 
# at Charité – Universitätsmedizin Berlin, with funding from the 
# DFG (FOR5657).
#####################################################################
"""


##################################################################
# Input:    - User-selected binary mask of coral skeleton
#           - User-inputted image resolution
#
# Output:   - Masks of total volume, total pore volume,
#             open pore volume, and closed pore volume
#           - Whole volume and per-slice quantification of
#             coral volume, total volume, total pore volume, 
#             open pore volume, closed pore volume, and derived
#             ratios
#        * Saved at input file path
##################################################################




import sys
import tkinter as tk
from tkinter import filedialog
import numpy as np
import tifffile as tif
from skimage import measure, morphology
import os
import pandas as pd
import matplotlib.pyplot as plt
import tifffile
from tkinter import simpledialog
from datetime import datetime
import threading
import itertools
import time



def spinner_task(stop_event):
    spinner = itertools.cycle(['|', '/', '-', '\\'])
    while not stop_event.is_set():
        sys.stdout.write(next(spinner))
        sys.stdout.flush()
        time.sleep(0.1)
        sys.stdout.write('\b')
        



######################################################################################
# Inputs
######################################################################################

######## User defined input
root = tk.Tk()
root.withdraw()  # Hide the root window
input_path = filedialog.askopenfilename(title="Select input binary image")

######## User defined pixel size
#root = tk.Tk()
#root.withdraw()
reso_um = simpledialog.askinteger(title="Pixel size", prompt="What is pixel size, in um?")
reso = reso_um*0.001 # units of mm


######################################################################################
# Processing
######################################################################################

######## Initialize output paths
input_filename = os.path.splitext(os.path.basename(input_path))[0]
path = os.path.dirname(input_path) 
outpath_whole = path + '/' + input_filename + '_results_whole.xlsx'
outpath_slices = path + '/' + input_filename + '_results_persclice.xlsx'


######## Extract input data (cropped binary image)
starttime = datetime.now()
print(str(starttime.strftime('%Hh:%Mm:%Ss')) + 'Beginning analysis of ' + input_filename)
print('Results will be saved in ' + path)
print('##################################')
stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Loading input data.          ", end="", flush=True)
#########################################
im = tifffile.imread(input_path)
imbin_raw = (im > 0).astype(np.uint8)
imbin = np.rollaxis(imbin_raw, 0,2)  # reshape to slice cross-sectionally
dz, dx, dy = np.shape(imbin)
#########################################
stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()
print()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Input data loaded')
print('Input size (x, y, z): ' + str(dx) + ' x ' + str(dy) + ' x ' + str(dz) + ' pixels')



######## Segmentation - pores and coral
stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Creating masks of pores, dissipments, and total volume.          ", end="", flush=True)
#########################################
# Pore and total (filled coral) masks
imbin_inv = 1-imbin # invert coral mask to get pores + background
#plt.imshow(imbin_inv[int(dz/2), :, :])
pores_and_bg_label = measure.label(imbin_inv) # label objects (pores + background)
#plt.imshow(pores_and_bg_label[int(dz/2), :, :])
unique, counts = np.unique(pores_and_bg_label[-1,:,:], return_counts=True) # get size (counts) of each object
ind = np.argmax(counts) # Identify label of largest object (which is background)
bg_label = unique[ind]
closed_pores = ~((pores_and_bg_label == 0) | (pores_and_bg_label == bg_label)) # remove background (including open pores) to get closed pores
#plt.imshow(closed_pores[int(dz/2), :, :])
ball20 = morphology.ball(radius=20)
total = morphology.binary_closing(imbin, ball20) # total volume via closing of coral mask
#plt.imshow(total[int(dz/2),:,:])
all_pores = np.bitwise_xor(total, imbin)
plt.imshow(all_pores[int(dz/2),:,:])
open_pores = np.bitwise_xor(all_pores, closed_pores)
#plt.imshow(open_pores[int(dz/2),:,:])
#########################################
stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Masks created')


######## Export created masks
stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Saving masks.          ", end="", flush=True)
#########################################
tif.imwrite(path + '/' + input_filename + '_closedpores.tiff', closed_pores.astype(np.int16))
tif.imwrite(path + '/' + input_filename + '_allpores.tiff', all_pores.astype(np.int16))
tif.imwrite(path + '/' + input_filename + '_openpores.tiff', open_pores.astype(np.int16))
tif.imwrite(path + '/' + input_filename + '_totalvolume.tiff', total.astype(np.int16))
#########################################
stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Masks saved')


######## Save representative slices of masks
fig = plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
im = plt.imshow(imbin[int(dz/2), :, :], cmap='Greys_r')
plt.title('Binary mask, coral, slice ' + str(int(dz/2)))
plt.tight_layout()
plt.savefig(path + '/' + input_filename + '_coral_mask_slc' + str(int(dz/2)) + '.png')
plt.close(fig)

fig = plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
im = plt.imshow(all_pores[int(dz/2), :, :], cmap='Greys_r')
plt.title('Binary mask, all pores, slice ' + str(int(dz/2)))
plt.tight_layout()
plt.savefig(path + '/' + input_filename + '_all_pores_mask_slc' + str(int(dz/2)) + '.png')
plt.close(fig)

fig = plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
im = plt.imshow(closed_pores[int(dz/2), :, :], cmap='Greys_r')
plt.title('Binary mask, closed pores, slice ' + str(int(dz/2)))
plt.tight_layout()
plt.savefig(path + '/' + input_filename + '_closed_pores_mask_slc' + str(int(dz/2)) + '.png')
plt.close(fig)

fig = plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
im = plt.imshow(open_pores[int(dz/2), :, :], cmap='Greys_r')
plt.title('Binary mask, open pores, slice ' + str(int(dz/2)))
plt.tight_layout()
plt.savefig(path + '/' + input_filename + '_open_pores_mask_slc' + str(int(dz/2)) + '.png')
plt.close(fig)

fig = plt.figure(figsize=(5,5))
plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
im = plt.imshow(total[int(dz/2), :, :], cmap='Greys')
plt.title('Binary mask, total volume, slice ' + str(int(dz/2)))
plt.tight_layout()
plt.savefig(path + '/' + input_filename + '_total_volume_mask_slc' + str(int(dz/2)) + '.png')
plt.close(fig)





######################################################################################
# Quantification
######################################################################################

######## Quantifying volumes
stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Quantifying vols.          ", end="", flush=True)
#########################################
# Whole volume quantification
vx = reso*reso*reso
whole_all_pore_vol = np.round(np.sum(vx*all_pores),5)
whole_closed_pore_vol = np.round(np.sum(vx*closed_pores),5)
whole_open_pore_vol = np.round(np.sum(vx*open_pores),5)
whole_coral_vol = np.round(np.sum(vx*imbin),5)
whole_total_vol = np.round(np.sum(vx*total),5)
whole_coral_to_total_perc = 100*whole_coral_vol/whole_total_vol
whole_closed_pore_perc = 100*whole_closed_pore_vol/whole_all_pore_vol
whole_total_porosity_perc = 100*whole_all_pore_vol/whole_total_vol
whole_closed_porosity_perc = 100*whole_closed_pore_vol/whole_total_vol
whole_open_porosity_perc = 100*whole_open_pore_vol/whole_total_vol

# Slice-wise volume quantification
slc_all_pore_vol = []
slc_closed_pore_vol = []
slc_open_pore_vol = []
slc_coral_vol = []
slc_total_vol = []
slc_coral_to_total_perc = []
slc_closed_pore_perc = []
slc_total_porosity_perc = []
slc_closed_porosity_perc = []
slc_open_porosity_perc = []
for i in range(dz):
    all_pv = np.round(np.sum(vx*all_pores[i,:,:]),5)
    closed_pv = np.round(np.sum(vx*closed_pores[i,:,:]),5)
    open_pv = np.round(np.sum(vx*open_pores[i,:,:]),5)
    cv = np.round(np.sum(vx*imbin[i,:,:]),5)
    tv = np.round(np.sum(vx*total[i,:,:]),5)
    
    slc_all_pore_vol.append(all_pv)
    slc_closed_pore_vol.append(closed_pv)
    slc_open_pore_vol.append(open_pv)
    
    slc_coral_vol.append(cv)
    slc_total_vol.append(tv)

    slc_coral_to_total_perc.append(np.round(100*cv/tv,2))
    slc_closed_pore_perc.append(np.round(100*closed_pv/all_pv,2))
    slc_total_porosity_perc.append(np.round(100*all_pv/tv,2))
    slc_closed_porosity_perc.append(np.round(100*closed_pv/tv,2))
    # open porosity = total porosity - closed porosity

#########################################
stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Outcomes quantified.')


stop_spinner = threading.Event()
spinner_thread = threading.Thread(target=spinner_task, args=(stop_spinner,))
spinner_thread.start()
print("Saving outcomes.          ", end="", flush=True)
#########################################
# Save vol results
df_vols = pd.DataFrame({'coral_vol_mm3': [whole_coral_vol],
                        'total_pore_vol_mm3': [whole_all_pore_vol], 
                        'closed_pore_vol_mm3': [whole_closed_pore_vol], 
                        'open_pore_vol_mm3': [whole_open_pore_vol],
                        'total_vol_mm3': [whole_total_vol], 
                        'coral_to_total_vol_perc': [whole_coral_to_total_perc],
                        'total_porosity_perc': [whole_total_porosity_perc],
                        'closed_porosity_perc': [whole_closed_porosity_perc],
                        'closed_to_all_pores_perc': [whole_closed_pore_perc]})
df_vols.T.to_excel(outpath_whole)
df_slc = pd.DataFrame({'coral_vol_mm3': slc_coral_vol[::-1],
                       'total_pore_vol_mm3': slc_all_pore_vol[::-1], 
                       'closed_pore_vol_mm3': slc_closed_pore_vol[::-1], 
                       'open_pore_vol_mm3': slc_open_pore_vol[::-1], 
                       'total_vol_mm3': slc_total_vol[::-1],
                       'coral_to_total_vol_perc': slc_coral_to_total_perc[::-1],                       
                       'total_porosity_perc': slc_total_porosity_perc[::-1],
                       'closed_porosity_perc': slc_closed_porosity_perc[::-1],
                       'closed_to_all_pores_perc': slc_closed_pore_perc[::-1]})
df_slc.index = range(1, dz+1)[::-1]
df_slc.to_excel(outpath_slices)   
#########################################
stop_spinner.set()
spinner_thread.join()
ctime = datetime.now()
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Outcomes saved')



endtime = datetime.now()
print('##################################')
print(str(ctime.strftime('%Hh:%Mm:%Ss')) + ': Analysis complete!')
duration = starttime - endtime
total_seconds = int(duration.total_seconds())
hours, remainder = divmod(total_seconds, 3600)
hours = hours + 1
minutes, seconds = divmod(remainder, 60)
print(f'Analysis duration = {hours}h:{minutes}m:{seconds}s')




