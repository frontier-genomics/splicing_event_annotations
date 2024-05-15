from src.eventAnnotateProcessor import EventAnnotate
from src.eventAnnotateProcessorList import EventAnnotateList

def run_workflow(input, dataset):

    if isinstance(input, list):  # If the input is a list, use BatchEventAnnotate
        processor = EventAnnotateList(input, dataset)
    else:  # Otherwise, use EventAnnotate
        processor = EventAnnotate(input['chrom'], input['start'], input['end'], input['strand'], input['transcript'], input['type'])

    return processor.process()