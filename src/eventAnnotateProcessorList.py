import logging
from src.eventAnnotateProcessor import EventAnnotate

logging.basicConfig(level=logging.INFO)

class EventAnnotateList:

    def __init__(self, inputs, dataset):
        self.inputs = inputs
        self.annotations = self.load_annotations(dataset)

    def load_annotations(self, dataset):
        first_input = self.inputs[0]
        print(first_input['chrom'], first_input['start'], first_input['end'], first_input['strand'], first_input['transcript'], first_input['type'])
        event_annotate = EventAnnotate(first_input['chrom'], first_input['start'], first_input['end'], first_input['strand'], first_input['transcript'], first_input['type'])
        print(event_annotate)
        self.annotations = event_annotate.read_refgene(dataset)
        print(self.annotations[0])

        return self.annotations
    
    def process(self):
        self.outputs = []
        for input in self.inputs:
            event_annotate = EventAnnotate(input['chrom'], input['start'], input['end'], input['strand'], input['transcript'], input['type'])
            output = event_annotate.process('refseq', get_annotations=self.annotations)
            self.outputs.append(output)
        return self.outputs