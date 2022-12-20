# SARS-CoV-2_spike_charge
SARS-CoV-2 spike charge analysis scripts

The repository contains data and scripts used in the manuscript "Evolution to increased positive charge on the viral spike protein may be part of  the adaptation of SARS-CoV-2 to human transmission"  by Cotten and Phan, currently under review at iScience and submitted to BioRxiv.


The following files are included:

Figure 1

Figure 1 data: Figure1_data.csv

Figure 1 Jupyter notebook: Figure1_VOCs_spike_first300_PCA_set4.ipynb

Figure 2

Figure 2 data can be found on Dropbox:

https://www.dropbox.com/s/pubmo5tu1rbznla/Complete_to_15Nov22_spike_charge_sample_total_charge_date.csv.zip?dl=0

The data file (too large to store on GitHub) can be retrieved from Dropbox, unzipped and loaded in the Figure2 Jupyter notebook to generate Figure 2.

Figure 2 Jupyter notebook: Figure_2_Complete_to_midNov22-spike_charge_for_GitHub.ipynb


Figure 3 

Data: 
cons_spike_from_first300_aln.fas (alignment of consensus spike sequences from first 300 genomes from select SARS-CoV-2 lineages)

SARS-CoV-2_SPIKE_features3.csv (key features of spike protein needed to make protein feature graph)

annotations_echt.csv (annotation table of protein substitutions)

Figure 3 script: HiLiter_spike_charge_GH.py

Command to run script: 

python HiLiter_spike_charge_GH.py cons_spike_from_first300_aln.fas SARS-CoV-2_SPIKE_features2.csv annotations_echt.csv

Figure 4, Panels A,B,C

Data:

229E_91_2022-1-27_spike_charge_table.csv

BCoV_114_114_2022-3-11_spike_charge_table.csv

OC43_235_2022-1-27_spike_charge_table.csv

PHEV_14_2022-3-26_spike_charge_table.csv

MERS-CoV_698_2022-1-27_spike_charge_table.csv

Figure 4, A,B,C Jupyter notebook: Figure4ABC_229E_OC43_MERS_Spike_charge.ipynb

Figure 4, Panels D,E, F

Data: 

Commands to run script:
Panel D

python HiLiter_spike_charge_GH_small.py 229E_human_camel_spike_con_aln.fas 229E_SPIKE_features.csv annotations_empty.csv

Panel E

python HiLiter_spike_charge_GH_small.py BCoV_OC43_cons_aln.fas OC43_SPIKE_features.csv annotations_empty.csv

Panel F

python HiLiter_spike_charge_GH_small.py PHEV_OC43_cons_aln.fas OC43_SPIKE_features.csv annotations_empty.csv





For questions or comments, contact 

Matthew Cotten Ph.D.
MRC Investigator
Professor of Pathogen Genomics and Bioinformatics
MRC/UVRI & London School of Hygiene and Tropical Medicine, Entebbe, Uganda 
Email: matthew.cotten@lshtm.ac.uk
