from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate
from splicing_event_annotator.eventAnnotateProcessorList import EventAnnotateList
from splicing_event_annotator.tsvAnnotateProcessor import TsvAnnotate
import csv

def run_workflow(input, dataset, genome, tsv = False, columns = [0,1,2,3,4,5], output_file = ""):

    if tsv == True:
        tsv_processor = TsvAnnotate(input, dataset, columns)
        input = tsv_processor.tsv

    if isinstance(input, list):
        processor = EventAnnotateList(input, dataset, genome)
    else:
        processor = EventAnnotate(input['chrom'], input['start'], input['end'], input['strand'], input['transcript'], input['type'])

    annotations = processor.process()

    if output_file != "":
        write_output(annotations, output_file)

    return annotations


def write_output(output, output_file):
    keys = output[0].keys()

    with open(output_file, 'w', newline='') as file:
        dict_writer = csv.DictWriter(file, keys, delimiter='\t')
        dict_writer.writeheader()
        dict_writer.writerows(output)