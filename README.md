# Downtime Assessment of Seismically Isolated Buildings

This page is the online repository for the following submitted manuscript: Shoma Kitayama, Enrique Abel Morales Moncayo, Anastasia Athanasiou. (Under Review). "Downtime Assessment of Seismically Isolated Buildings Designed by the Procedures of ASCE/SEI 7-16."
 
This repository contains downtime assessment files programmed in MATLAB (version 9.10.0.1649659 (R2021a) Update 1; The MathWorks, Inc).

Prepared by:
Shoma Kitayama (s.kitayama@leeds.ac.uk) - University of Leeds

Checked by:
Enrique Abel Morales Moncayo (eamorales5@espe.edu.ec) - Universidad de las Fuerzas Armadas (ESPE); Anastasia Athanasiou (anaj.athanasiou@gmail.com) - Concordia University

# Log

27.Mar.2022, Version: 2 (Last modified)

04.Oct.2021, Version: 1

# Notes

The folder contains MATLAB codes that computes the following values that are the results of seismic downtime assessment of seismically isolated buildings with SCBF designed by RI=2 and TFP-1 based on Conditional Spectra approach (notations are explained in the manuscript):

- Downtime vulnerability curves
- Mean mobilization times
- Expected annual downtime (EADT)
- Expected annual downtime losses (EADTL)
 
To implement this example code, you need to download the data of results of nonlinear seismic response analysis from the following link:

https://drive.google.com/file/d/1pNAoIicndCVjONx4P9fMqTAlJawJbsfy/view?usp=sharing

The downloaded data (i.e., a folder, after extraction of zip file, "Response_analysis_data_CS_SCBF_RI2.0_Iso1_ASCE16") of results of seismic respone analysis should be located to the folder in the same place as the matlab files are located.

The zip file should be extracted to get a folder that contains data. The data is the results of analysis of seimically isolated building with SCBF (superstructure is designed per RI=2.0) and with the smallest required isolators (TFP-1) based on ASCE7-16 (2017) seismic design standard.

Please note that one of the MATLAB file "fn_mle_pc.m" used in this folder was obtained from the data set that is associated with the following journal article: Baker JW. (2015). “Efficient analytical fragility function fitting using dynamic structural analysis.” Earthquake Spectra, 31(1), 579-599.
