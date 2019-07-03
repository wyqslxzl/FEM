# FEM
paper code 
“Multi-Region Brain Stimulation Optimization Using Transcranial Direct Current Stimulation”

# m files
DC_FE_ME_S:  Brain stimulation optimization using ME scheme in single target case.
DC_FE_ME_M:  Brain stimulation optimization using weighted ME scheme in multi-target case.
DC_FE_MLS_S:  Brain stimulation optimization using MLS scheme in single target case.
DC_FE_MLS_M:  Brain stimulation optimization using weighted MLS scheme in multi-target case.

fun_ME_S: cost funtion of ME scheme.
fun_ME_M: cost funtion of weighted ME scheme.
fun_MLS_S: cost funtion of MLS scheme.
fun_MLS_M: cost funtion of weighted MLS scheme.

con_ME_S:  constrains of ME scheme.
con_ME_M:  constrains of weighted ME scheme.
con_MLS_S:  constrains of MLS scheme.
con_MLS_M:  constrains of weighted MLS scheme.

# HC001.rar
The files related to head model, including HC001_AAL.nii, HC001_Atlas.nii, HC001_brain_seg_1.nii and HC001_brain_seg_2.nii.
HC001_AAL.nii: AAL altas mapped gray matter template.
HC001_Atlas.nii: SWMT task T-map mapped gray matter template.
HC001_brain_seg_1.nii: gray matter template.
HC001_brain_seg_2.nii: white matter template.

# DC_64 folder
DC_64 folder contains DC_64.mat file which saves the coeffiecint matrix A of 63 freedom electrode.Because this file has 1.24G, so we compress it into 64 20M parts.
