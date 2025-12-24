
import sys
from Bio import SeqIO
import re
import argparse
import pandas as pd
import os

def parse_inputs(argument_list):
	"""
	Function: 
		This function processes the command-line arguments passed to the script using sys.argv[]

	Parameters: 
		argument_list: command-line arguments, uses sys.argv[1:] (excludes the script name at sys.argv[0])
	
	Returns:
		args: parsed arguments list
	"""
	parser = argparse.ArgumentParser(description="Generates a csv file with the counts of occurences of primer template combinations in sequences, for a specified experiment.")

	parser.add_argument('--run', required=True, help='CSV file with all the SRR run numbers:')
	parser.add_argument('--meta', required=True, help='CSV file with the metadata table (where experiment information such as tempature is)')
	parser.add_argument('--exp', required=True, type=str, help='Name of the experiment of interest (example: B27)')
	parser.add_argument('--map', required=True, help='CSV file with mapping information for generation of the count file')
	parser.add_argument('--t', required=True, help='A txt file with the new template names(ie: ST0, ST1) of your desired experiment. Delinated by a space.')
	parser.add_argument('--p', required=True, help='A txt file with the new primer names(ie: V0, V1) of your desired experiment. Delinated by a space.')

	args = parser.parse_args(argument_list)
	
	return args


def get_metadata(run_table, metadata, exp_name):
	"""
	Function: 
		This function creates a pandas data frame with the metadata (temperature, PCR method etc.) and SRR run IDs for the experiment of interest.

	Parameters:
		run_table (.csv): csv file with the SRR run IDs.
		metadata (.csv): csv file with all metadata.
		exp_name (str): the experiment name of interest (for example: 'B27' or 'A10')

	Returns:
		filtered_meta (pandas dataframe): with the metadata and SRR run IDs for the specified experiment.
	""" 	
	run_table = pd.read_csv(run_table)
	metadata = pd.read_csv(metadata)

	table = run_table[['Run', 'BioSample']]
	
	metadata = metadata.rename(columns={'Biosample Accession': 'BioSample'})
	exp_metadata = metadata[metadata['Current Exp Name'].str.contains(exp_name)]

	filtered_meta = pd.merge(table, exp_metadata, on='BioSample')
	
	return filtered_meta.reset_index(drop=True)

def get_run_IDs(metadata):
	"""
	Function: 
		This function creates a list of SRR run IDs with the file extension ".fastq", that correspond to the desired experiment. This list is used when extracting sequences from the SRR files.

	Parameters: 
		metadata (pandas dataframe): output from get_metadata.

	Returns: 
		run_IDs_list (list): list of all SRR run IDs with the sting ".fastq" added at the end of each ID.
	"""
	temp_list = list(metadata['Run'])
	file_type = ".fastq"
	run_IDs_list = [item + file_type for item in temp_list]
	
	return run_IDs_list

def get_templates(templates_file):
	"""
	Function: 
		This function extracts the template names from a .txt file for the specified experiment.

	Parameters: 
		templates_file (.txt): User created .txt file with SPACE DELINATED template names. Reference "Table 1" from the study to find which templates correspond to your experiemnt of interest. 

	Example .txt file contents: 
		ST0 ST1 ST3 ST4 ST5 ST6 ST7 ST8 ST9

	Returns: 
		templates_list (list): List of template names
	"""
	with open(templates_file, 'r') as file:
		templates = file.read().split()
	
	return templates

def get_primers(primers_file):
	"""
	Function: 
		This function extracts the primer names from a .txt file for the specified experiment.

	Parameters: 
		primers_file (.txt): User created .txt file with SPACE DELINATED primer names. Reference "Table 1" from the study to find which primer correspond to your experiemnt of interest. 

	Example .txt file contents: 
		V10 V11 V12 V13 ... V36

	Returns: 
		primers_list (list): List of primer names
	"""
	with open (primers_file, 'r') as file:
		primers = file.read().split()
	
	return primers

def get_mapping_info(map, templates, primers):
	"""
	Function: This function filters the rows from the mapping file using the templates_list and primers_list from the functions get_templates and get_primers, to obtain only the rows that correspond to the experiemnt of interest. 

	Parameters: 
		mapFile (.csv): mapping file from supplementary data in the study (TableS3.xlsx) Note: you must change the file type from .xlsx to .csv.
		templates (list): output list from the function get_templates.
		primers (list): output list from the function get_primers.

	Returns: 
		filtered_map (pandas dataframe): filtered to only contain rows corresponding to the templates and primers listed in the templates and primers lists.
	"""
	map = pd.read_csv(map, sep=None, engine='python')
	filtered_map = map[map['current_template_name'].isin(templates) & map['current_primer_name'].isin(primers)]
	filtered_map = filtered_map.reset_index(drop=True)

	return filtered_map

def get_sequences(run_IDs_list):
	"""
	Function:
		The function uses the list run_IDs_list, which is the output from the function get_run_IDs and extracts all sequences from the specified SRR run IDs in the run_IDs_list.  
		Note: The file path to where the sequencing files are is hardcoded. The codes is designed to run on from the users directory within the server. 

	Parameters:
		run_IDs_list (list): output list from the function get_run_IDs, which contains all SRR run IDs for the desired experiement (because the extension .fastq has been added, this list is the list of file names)

	Returns: 
		seq_dict (dictionary): key{SRR run IDs}, value{sequences}
			Example: key{SRR8399749.fastq}, value{ATCTTGCCGTAAGTC....}
	"""
	filepath = '/home/mpa' #hardcoded filepath 
	seq_dict = {}
	for file in os.listdir(filepath):
		if file in run_IDs_list:
			full_path = os.path.join(filepath, file)
			seq_dict[file] = []

			with open(full_path) as handle:
				for record in SeqIO.parse(handle, "fastq"):
					seq_dict[file].append(str(record.seq))			
	
	return seq_dict			
 
def get_counts(filtered_map, seq_dict):
	"""
	Function: 
		This function counts the number of occurences of each primer template combination found in all the sequences from seq_dict.
		Note: column names are hardcoded

	Parameters:
		filtered_map (pandas dataframe): output from get_mapping_info, which contains the rows of template and primer combinations, filtered to only contains those combinations possible for the specified experiment. 
		seq_dict (dictionary): output from the function get_sequences which contains the following keys and values: key{SRR run IDs}, value{sequences}

	Returns: 
		counts_pivot (pandas dataframe): pivot table with the "feature_name" as the index, the 'run_ID' as the columns and the values as the "counts", where the counts are the summed number of occurances of each primer and template combination found that SRR file.
			
			Example output (numbers inaccurate, just an example): 
			____________________________________________________________________________________
				R_feature_name	SRR8399749.fastq  SRR83399769.fastq .......... SRR8399829.fastq
					ST01V10				12				70								500
					ST01V11			    2				100								10
					.					.				.								.
					.					.				.								.
					.					.				.								.
					ST09V36				0				25								45
			____________________________________________________________________________________
	"""
	table = []
	for i, row in filtered_map.iterrows():
		template_pattern = row['template_recognition_sequence']
		primer_pattern = row["Primer Seq (5'-3')"]
		template_primer = row["R_feature_name"]

		for run_ID, seqs in seq_dict.items():
			for seq in seqs:
				if template_pattern in seq and primer_pattern in seq:
					table.append({
						'R_feature_name': template_primer,
						'run_ID': run_ID,
						'count': 1
					})
	
	counts_df = pd.DataFrame(table)
	counts_df = counts_df.groupby(['R_feature_name', 'run_ID'])['count'].sum().reset_index()

	counts_pivot = counts_df.pivot(
		index='R_feature_name', 
		columns='run_ID', 
		values='count'
		).fillna(0).astype(int)
	
	return counts_pivot

def combine_counts_meta(counts_pivot, filtered_meta): 
	"""
	Function: 
		This function combines the counts_pivot table result from the get_counts function with the filtered_meta from get_metadata function. 

	Parameters:
		counts_pivot (pandas dataframe): output table from the get_counts function, where the table contains the counts of occurences of primer template combinations in each sequence of a run file (SRR).
		filtered_meta (pandas dataframe): output table from the get_metadata function, where the result is a table with only the metadata of the desired experiment. 

	Returns: 
		counts_meta_table (pandas dataframe): table with counts and metadata.
		
			Example output table (not accurate info, just an example): 
				__________________________________________________________________________________
					Rpl		Run			R_feature_name 	Annealing Temp	PCR Method	Count
					1		SRR8399877		ST09V34					55C				TAS		50
					1		SRR8399795		ST09V30					45C				DePCR	100
					.			.				.					.				.		.
					.			.				.					.				.		.
					8		SRR8399793		ST09V33					55C				DePCR	0
				____________________________________________________________________________________
	"""
	counts_long = counts_pivot.reset_index().melt(id_vars='R_feature_name',
                               var_name='Run',
                               value_name='Count'
                              )
    
	counts_long['Run'] = counts_long['Run'].str.replace(r'\.[^.]+$', '', regex=True)
    
	filtered_meta = filtered_meta[['Run', 'PCR Method', 'Annealing Temp', 'Rpl']]
	combined_df = pd.merge(counts_long, filtered_meta, on='Run')
	new_order = ['Rpl', 'Run', 'R_feature_name', 'Annealing Temp', 'PCR Method', 'Count']
	counts_meta_table = combined_df[new_order]
	
	return counts_meta_table

def mismatch_matrix(counts_meta_table, filtered_map):
	"""
	Function:
		This function creates columns in the output table counts_meta_table from the function combine_counts_meta, titled 5', Middle, and 3' with values as 0 (no mismatch) or 1 (mismatch). These columns are a matrix that say whether there is a mismatch or no mismatch at the 5', Middle and 3' positions in the primer sequence.

	Parameters:
		counts_meta_table (pandas dataframe): output from the function combine_counts_meta function.
		filtered_map (pandas dataframe): output from the function filtered_map.

	Returns: 
		final_table (pandas dataframe): table with the counts, metadata, and mismatch matrix
			Example output table (not accurate info, just an example): 
				_________________________________________________________________________________________________
					Rpl		Run			R_feature_name			5'	Middle	3' 	Annealing Temp	PCR Method	Count
					1		SRR8399877		ST09V34				0		1	0		55C				TAS		50
					1		SRR8399795		ST09V30				1		0	1		45C				DePCR	100
					.			.				.				.		.	.		.				.		.	
					.			.				.				.		.	.		.				.		.
					8		SRR8399793		ST09V33				0		0	0		55C				DePCR	0
				__________________________________________________________________________________________________
		
		
	""" 
	perfect_match_data = {'templates': ['ST00', 'ST01', 'ST02', 'ST03', 'ST04', 'ST05', 'ST06', 'ST07', 'ST08', 'ST09'],
					   'perfect_match': ['CTA', 'GTA', 'TTA', 'ATA', 'CAA', 'CCA', 'CGA', 'CTC', 'CTG', 'CTT']}
	df = pd.DataFrame(perfect_match_data)
		
	temp_map = filtered_map[['R_feature_name', "variable positions -14,-8,-2 from 3' end"]]
	final_map = temp_map.rename(columns={"variable positions -14,-8,-2 from 3' end": 'variable'})
		
	df2 = pd.merge(counts_meta_table, final_map, on='R_feature_name')
		
	df2.insert(loc=3, column="5'", value=0)
	df2.insert(loc=4, column="Middle", value=0)
	df2.insert(loc=5, column="3'", value=0)

	for row in range (len(df2)):
		for row2 in range (len(df)):
			if df2.loc[row, 'R_feature_name'][:4] == df.loc[row2, 'templates'][:4]:
				str1 = df.loc[row2, 'perfect_match']
				str2 = df2.loc[row, 'variable']
				
				df2.loc[row, "5'"] = 1 if str2[0] != str1[0] else 0
				df2.loc[row, "Middle"] = 1 if str2[1] != str1[1] else 0
				df2.loc[row, "3'"] = 1 if str2[-1] != str1[-1] else 0
	
	df2 = df2.drop('variable', axis=1)
	final_table = df2.set_index('Rpl').sort_index()

	final_table.to_csv('final.csv', index=True)
	
	
	return final_table

def main():

    #Parses inputs from command-line.
	args = parse_inputs(sys.argv[1:])

	#Creates a table with the metadata for the specific experiment.
	filtered_meta = get_metadata(args.run, args.meta, args.exp)

	#Extracts SRR run ids from output of get_metadata and stores in a list.
	run_IDs_list = get_run_IDs(filtered_meta)

	#Extracts sequences and stores in a dictionary with the SRR run number as the key and the sequences as the value, using the list from get_run_IDs.
	seq_dict = get_sequences(run_IDs_list)

	#Gets template names from a user made text file and stores in a list.
	templates_list = get_templates(args.t)

	#Gets primer names from a user made text file and stores in a list.
	primers_list = get_primers(args.p)

	#Gets map_table with only desired templates and primers from resulting list from the get_templates and get_primers functions.
	filtered_map = get_mapping_info(args.map, templates_list, primers_list)

	#Calculates the counts of each primer template combo from the get_mapping_info function in the sequences from the get_sequences function.
	counts_pivot = get_counts(filtered_map, seq_dict)

	#Combines the filtered_meta table from the get_metadata function and the counts_pivot table from the get_counts function
	counts_meta_table = combine_counts_meta(counts_pivot, filtered_meta)

	#adds the mismatch matrix to the output table from counts_meta_table and saves to a .csv file titled "final.csv"
	final_table = mismatch_matrix(counts_meta_table, filtered_map)

if __name__ == "__main__":
	main()
