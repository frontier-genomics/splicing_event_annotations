from src.eventAnnotateProcessor import EventAnnotate
from src.eventAnnotateProcessorList import EventAnnotateList
from src.tsvAnnotateProcessor import TsvAnnotate
import csv

def run_workflow(input, dataset, tsv = False, columns = [0,1,2,3,4,5]):

    if tsv == True:
        tsv_processor = TsvAnnotate(input, dataset, columns)
        input = tsv_processor.tsv

    if isinstance(input, list):  # If the input is a list, use BatchEventAnnotate
        processor = EventAnnotateList(input, dataset)
    else:  # Otherwise, use EventAnnotate
        processor = EventAnnotate(input['chrom'], input['start'], input['end'], input['strand'], input['transcript'], input['type'])

    return processor.process()

def write_output(output, output_file):
    # Get the keys (column names) from the first dictionary in the list
    keys = output[0].keys()

    # Write list of dictionaries to CSV
    with open(output_file, 'w', newline='') as file:
        dict_writer = csv.DictWriter(file, keys, delimiter='\t')  # Set delimiter to '\t' for TSV
        dict_writer.writeheader()
        dict_writer.writerows(output)