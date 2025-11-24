import logging
from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate

class EventAnnotateList:

    def __init__(self, inputs, dataset, genome):
        self.inputs = inputs
        self.dataset = dataset
        self.genome = genome
        self.annotations = self.load_annotations(dataset, genome)

    def load_annotations(self, dataset, genome):
        first_input = self.inputs[0]
        event_annotate = EventAnnotate(first_input['chrom'], first_input['start'], first_input['end'], first_input['strand'], first_input['transcript'], first_input['type'])
        self.annotations = event_annotate.read_refgene(dataset, genome)
        logging.info(f"{dataset} annotations loaded")

        return self.annotations
    
    def process(self):
        logging.info(f'Processing {len(self.inputs)} inputs')
        self.outputs = []
        
        for i, input in enumerate(self.inputs):
            index = i + 1
            
            # Log every 1000 records
            if index % 1000 == 0:
                logging.info(f'Processed {index} records')

            event_annotate = EventAnnotate(input['chrom'], input['start'], input['end'], input['strand'], input['transcript'], input['type'])
            output = event_annotate.process(self.dataset, self.genome, get_annotations=self.annotations)
            self.outputs.append(output)

        logging.info(f'Finished processing {len(self.inputs)} inputs')

        return self.outputs