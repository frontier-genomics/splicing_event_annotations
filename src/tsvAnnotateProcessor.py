import pandas as pd
import numpy as np
import src.eventAnnotateProcessor as eap

class TsvAnnotate:
    def __init__(self, tsv_path, input_columns = [0,1,2,3], tx_column = 999, sep='\t'):
        self.tsv = pd.read_csv(tsv_path, sep=sep)
        print(f"Opening tsv file: {tsv_path}")
        print(f"Rows: {self.tsv.shape[0]}")
        event_number = range(self.tsv.shape[0])
        print(event_number)
        print(f"Found {max(event_number)+1} events in tsv file")
        chrom = self.tsv.iloc[:, input_columns[0]]
        start = self.tsv.iloc[:, input_columns[1]]
        end = self.tsv.iloc[:, input_columns[2]]
        strand = self.tsv.iloc[:, input_columns[3]]
        transcript = self.tsv.iloc[:, tx_column] if tx_column != 999 else ["NA"] * (max(event_number) + 1)

        self.events = {'chrom': chrom,
                       'start': start,
                       'end': end,
                       'strand': strand,
                       'transcript': transcript,
                       'type': ["sj"] * (max(event_number) + 1),
                       'event_number': event_number}

    def annotate(self, dataset = 'refseq'):
        self.annotations = []

        for index in self.events['event_number']:

            print(f"Annotating event {index}...")

            event = eap.EventAnnotate(chrom = self.events['chrom'][index],
                                      start = self.events['start'][index],
                                      end = self.events['end'][index],
                                      strand = self.events['strand'][index],
                                      transcript = self.events['transcript'][index],
                                      type = self.events['type'][index])
            
            # if event.coordinates['transcript'] == "get transcripts":
            #     transcript = event.get_mane_transcript(annotation)
            
            print(dataset)

            event.get_annotations(dataset)

            start = event.reference_match('start')
            end = event.reference_match('end')

            annotation = event.fetch_transcript_annotations(start, end)['event']

            print(annotation)

            self.annotations.append(annotation)

        return self.annotations