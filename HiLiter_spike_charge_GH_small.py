#!/usr/local/bin/python
import sys
import os.path
from Bio import SeqIO
import csv
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict
 
if len(sys.argv)!= 4:
	print("Usage: python HiLiter_spike_charge_GH.py sequences.fasta feature_table.csv annotation_table.csv")
	sys.exit()
	
print("Running HiLiter_spike_charge_GH.py really fast")
all_seqs= sys.argv[1]
outprefix = os.path.splitext(all_seqs)[0]
feature_table = sys.argv[2]
annotation_table = sys.argv[3]
protein_name = "Spike"
annotate_cutoff = 1
seq_names =[]
sequences =[]

for record in SeqIO.parse(open(all_seqs, "rU"), "fasta"):
	seq_names.append(record.id)
	sequence_string = str(record.seq)	
	sequences.append(sequence_string)
total_number_sequences = len(sequences)
genome_length = len(sequences[0])
ytick_label_size = 5
genome_length = len(sequences[0])

marker_size = 6	
graph_pad = 10

graph_names=[]
graph_names.append(seq_names[0])
neg_charge_list = ['D', 'E']
pos_charge_list = ['K','R']
neutral_charge_list = ['A','C','F','G','H','I','L','M','N','P','Q','S','T','V','W','Y']

#Function to generate lists of differences between ref and test sequences
def vergleichen(seqA, seqB, nameB, index):	
	for i in range (len(seqA)):
		if seqA[i] in neg_charge_list and seqB[i] in neutral_charge_list:
			color = "#ff5c33" #orange
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] in neg_charge_list  and seqB[i] in pos_charge_list:
			color = "crimson" 
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] in neutral_charge_list and seqB[i] in neg_charge_list:
			color = "#538cc6" #less dk blue
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] in neutral_charge_list and seqB[i] in pos_charge_list:
			color = "#ff5c33" #orange
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] in pos_charge_list and seqB[i] in neutral_charge_list:
			color = "#538cc6" #less dk blue
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] in pos_charge_list and seqB[i] in neg_charge_list:
			color = "#204060" #dk blue
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] in neg_charge_list and seqB[i] =="-":
			color = "#ff5c33" #orange
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] in pos_charge_list and seqB[i] =="-":
			color = "#538cc6" #less dk blue
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] == "-" and seqB[i] in neg_charge_list:
			color = "#538cc6" #less dk blue
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)
		elif seqA[i] == "-" and seqB[i] in pos_charge_list:
			color = "#ff5c33" #orange
			line = index
			diff_list.append((i+1,line,color))
			position = str(i+1)
			change_annotation = str(seqA[i]+position+seqB[i])
			table_list.append((index,nameB,i+1,seqA[i],seqB[i],change_annotation))
			complete_change_list.append(i+1)			
		else:
			continue
diff_list=[]
table_list=[]
complete_change_list = []
 
#call vergleichen on each pair
for i in range(1, len(sequences)):
	vergleichen(sequences[0], sequences[i], seq_names[i], i)

#generate table with all changes, format in format to be used by annotation table
all_changes=[]
difference_summary = open(outprefix +"_all_changes.csv", "a")
for element in table_list:
	row_number =  str(element[0])
	sample = element[1]
	position = str(element[2])
	change_annotation = str(element[3]+position+element[4])
	print_string = str(sample+","+row_number+","+position+","+change_annotation)
	print(print_string, file=difference_summary)
	all_changes.append(change_annotation)
difference_summary.close()

#generate table with all changes for single bar graph
occurrence = defaultdict(lambda: 0)

for element in all_changes:
	occurrence[element] += 1

#generate annotation list
annotation_list=[]
frequency_count = defaultdict(int)
for item in all_changes:
	frequency_count[item] += 1
for key, value in list(frequency_count.items()):
	if "-" in key:
		pass
	elif value >= annotate_cutoff:
		annotation_list.append(key)

def getKey(item):
	return item[0]
	
#now generate graph of differences. 
ticks=[] 
for x in range(len(seq_names)):
	ticks.append(x)
		
fig = plt.figure()

#ax2 = plt.subplot2grid((28,1), (6,0), rowspan = 22)
ax2 = plt.subplot2grid((28,1), (6,0), rowspan = 7)

ax2.set_ylim(0.9, len(sequences)) 

ax2.set_xlim(-graph_pad,(genome_length+graph_pad))

ax2.set_yticks(ticks) 
ax2.set_yticklabels(seq_names, size = ytick_label_size)
ax2.grid(False)
ax2.set_xlabel(protein_name+' Protein Position')

for i in range(len(diff_list)):
	color = diff_list[i][2]
	ax2.broken_barh([(int(diff_list[i][0]), marker_size)] , (int(diff_list[i][1]), 0.85), facecolors=(color),edgecolors=('None')) 

#annotate the important aa changes. 
annotation_table_rdr = csv.reader(open(annotation_table, "rU"))
for row in annotation_table_rdr:
	x_pos = int(row[2])+7
	y_pos = float(row[0])
	change = str(row[5])
	ax2.annotate(change, xy=(1,1), xytext=(x_pos, y_pos), fontsize=ytick_label_size)	

#add genome feature graph
feature_table = csv.reader(open(feature_table, "rU"))
feature_names =[]
feature_map_positions = []
feature_lengths = []
feature_colors=[]
feature_row = []
next(feature_table)#skip header
for row in feature_table:	
	feature_name = str(row[0])
	feature_names.append(feature_name)
	feature_start = int(row[1])
	feature_map_positions.append(feature_start)
	feature_length = int(row[2])- int(row[1])
	feature_lengths.append(feature_length)
	feature_color = str(row[4])
	feature_colors.append(feature_color)
	feature_row.append(int(row[3]))
	
ax3 = plt.subplot2grid((28,1), (1,0), rowspan = 5)

ax3.set_ylim(0.8, 6.2) 
ax3.set_xlim(-graph_pad,(genome_length+graph_pad))
ax3.set_xticklabels([])
ax3.set_yticklabels([])

ax3.grid(False)
ax3.get_yaxis().set_visible(False)

for i in range(len(feature_map_positions)):
	ax3.broken_barh([(int(feature_map_positions[i]), int(feature_lengths[i]))],(float(feature_row[i]), 0.9), facecolors=(str(feature_colors[i])), edgecolors=('None'))
	ax3.annotate(feature_names[i], xy=(1,1), xytext=(feature_map_positions[i]+int(feature_lengths[i])+5, (float(feature_row[i]))+0.2), fontsize=5)

plt.savefig(outprefix+'_differences_new.pdf', bbox_inches='tight', dpi=300)
plt.savefig(outprefix+'_differences_new.jpg', bbox_inches='tight', dpi=300)

print("That's All Folks!")
