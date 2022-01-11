#!/usr/bin/env python

import os,glob,gzip,shutil
from fpdf import FPDF
from datetime import datetime

#Initialize pdf
pdf=FPDF()
pdf.add_page()
page_width = pdf.w - 2 * pdf.l_margin
pdf.set_font('Times','B',16.0)
pdf.cell(page_width,0.0,'SCZ White Matter Pipeline Correlation Matrix',align='C')
pdf.ln(8)

# Write xnat specific project info to header
pdf.set_font('Times','',14.0)
project = os.environ["xnat_project"]
subject = os.environ["xnat_subject"]
session = os.environ["xnat_session"]
fmri_scan = os.environ["fmri_niigz"]

scan = fmri_scan.split('/')[-1]

pdf.cell(page_width,0.0,'XNAT Project: ' + project + ' Subject: ' + subject + ' Session: ' +  session,align='C')
pdf.ln(8)

# Get current time and write to header
now = datetime.now()
dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

pdf.cell(page_width,0.0,dt_string,align='C')
pdf.ln(12)

#Load in png file and write to pdf
first = True
for png_file in glob.glob('result1_corrmatrix/matr_*.png'):
	
	# To do: Split scan name and add check for multiple fmri inputs


	#Don't add page if we are on the first image
	if first:
		first = False
	else:
		pdf.add_page()
	pdf.set_font('Times','B',14.0)
	pdf.cell(page_width,0.0,scan)
	pdf.ln(6)
	pdf.image(png_file,0,None,200,0)      
	#pdf.ln(20)
	#pdf.add_page()

# Finalize pdf
pdf.output(project + '_' + subject + '_SCZ_CorrMatr.pdf','F')

# Compress nifti files
for nii in glob.glob('**/*.nii',recursive=True):
	with open(nii,'rb') as f_in, gzip.open(nii + '.gz', 'wb') as f_out:
		f_out.writelines(f_in)
	os.unlink(nii)
