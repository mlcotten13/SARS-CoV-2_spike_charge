# SARS-CoV-2_spike_charge
SARS-CoV-2 spike charge analysis scripts

The repository contains data and scripts used in the manuscript "Evolution to increased positive charge on the viral spike protein may be part of  the adaptation of SARS-CoV-2 to human transmission"  by Cotten and Phan, currently under review at  Virus Evolution and submitted to BioRxiv.

The following files are be included:

Figure 1
Figure 1 data: Figure1_data.csv
Figure 1 Jupyter notebook: Figure1_VOCs_spike_first300_PCA_set4.ipynb

Figure 2
Figure 2 data can be found on Dropbox:
https://www.dropbox.com/s/ii1uf0fzyxem3h1/Figure2_data.csv.zip?dl=0
The data file (too large to store on GetHub) can be retrieved from Dropbox, unzipped and loaded in the Figure2 Jupyter notebook to generate Figure 2.
Figure 2 Jupyter notebook: Figure2_Complete_to_midJune22-spike_charge.ipynb

Figure 3 
Data: 
cons_spike_from_first300_aln.fas (alignment of consensus spike sequences from first 300 genomes from select SARS-CoV-2 lineages)
SARS-CoV-2_SPIKE_features3.csv (key features of spike protein needed to make protein feature graph)
annotations_echt.csv (annotation table of protein substitutions)
Figure 3 script: HiLiter_spike_charge_GH.py
command to run script: 
python HiLiter_spike_charge_GH.py cons_spike_from_first300_aln.fas SARS-CoV-2_SPIKE_features2.csv annotations_echt.csv

Figure 4
Data:

Figure 4 Jupyter notebook


For questions or comments, contact 

Matthew Cotten Ph.D.
MRC Investigator
Professor of Pathogen Genomics and Bioinformatics
MRC/UVRI & London School of Hygiene and Tropical Medicine, Entebbe, Uganda 
Email: matthew.cotten@lshtm.ac.uk
