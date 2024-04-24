import pandas as pd
import numpy as np

class EventAnnotate:

    def __init__(self, chrom, start, end, strand, type, annotation_choice):
        self.coordinates = {
            'chrom': str(chrom),
            'start': np.int64(start),
            'end': np.int64(end),
            'strand': str(strand),
            'type': str(type)
        }
        self.annotation_choice = annotation_choice

    def get_annotations(self):
        if self.annotation_choice == "refseq_mane":
            self.annotation = pd.read_csv("resources/annotations/refseq_mane_introns_sorted.tsv", sep='\t')
            self.exons = pd.read_csv("resources/annotations/refseq_mane_exons_sorted.tsv", sep='\t')
        elif self.annotation_choice == "refseq_curated":
            self.annotation = pd.read_csv("resources/annotations/refseq_curated_introns_sorted.tsv", sep='\t')
            self.exons = pd.read_csv("resources/annotations/refseq_curated_exons_sorted.tsv", sep='\t')
        else:
            raise ValueError("Invalid annotation choice. Please select either 'refseq_mane' or 'refseq_curated'.")
        
        self.annotation.columns = ['chrom', 'start', 'end', 'transcript', 'intron', 'strand']

    def matcher(self, start_end, base = 1):
        if start_end not in ['start', 'end']:
            raise ValueError("Invalid start_end choice. Please select either 'start' or 'end'.")
        
        if start_end == "start":
            query = self.coordinates[start_end] - base
        elif start_end == "end":
            query = self.coordinates[start_end]
        
        chrom = self.coordinates['chrom']
        
        matching_rows = self.annotation[
            (self.annotation[start_end] == query) &
            (self.annotation['chrom'] == chrom)
        ]
        
        matching_rows['start'] = matching_rows['start'] + base
        match = pd.DataFrame(matching_rows)
        return(match)
        
    def create_annotations(self, start_matches, end_matches):
        if end_matches.empty and start_matches.empty:
            return {'transcript': "unknown", 'event': "unannotated junctions"}
        
        elif start_matches.empty or end_matches.empty:
            if end_matches.empty:
                print("start matches")
                end_matches = ""
                start_intron = start_matches['intron'].iloc[0]
                start_tx = start_matches['transcript'].iloc[0]
                distance = self.coordinates['end'] - start_matches['end'].iloc[0]
                return(self._name_event(start_tx, distance, self.coordinates['end'], "acceptor", start_matches))
            
            elif start_matches.empty:
                print("end matches")
                start_matches = ""
                end_intron = end_matches['intron'].iloc[0]
                end_tx = end_matches['transcript'].iloc[0]
                distance = end_matches['start'].iloc[0] - self.coordinates['start']
                return(self._name_event(end_tx, distance, self.coordinates['start'], "donor", end_matches))
        
        else:
            start_intron = start_matches['intron'].iloc[0]
            end_intron = end_matches['intron'].iloc[0]
            start_tx = start_matches['transcript'].iloc[0]
            end_tx = start_matches['transcript'].iloc[0]

            if start_tx == end_tx:
                if self.coordinates['type'] == "ir":
                    event = f"intron {start_intron} retention"
                    return {'transcript': start_tx, 'event': event}
                elif start_intron == end_intron:
                    event = f"canonical exon {start_intron}-{end_intron+1} splicing"
                    return {'transcript': start_tx, 'event': event}
                elif start_intron != end_intron:
                    event = f"exon {'-'.join(str(i) for i in range(start_intron+1, end_intron+1))} skipping"
                    return {'transcript': start_tx, 'event': event}

            elif start_tx != end_tx:
                    return {'transcript': "transcript mismatch", 'event': "unknown event"}
                
    def _name_event(self, tx, distance, start_end, splice_site_type, matches):
        direction = "+" if distance > 0 else ""
        within_tx_intron = self.annotation[(self.annotation['transcript'] == tx) &
                        (self.annotation['start'] <= start_end) &
                        (self.annotation['end'] >= start_end)
        ]
        
        if not within_tx_intron.empty:
            print("found an intronic variant")
            print(f"{start_end} is between \n {within_tx_intron}")
            print(within_tx_intron['intron'].iloc[0])
            print(matches['intron'].iloc[0])
            
            if within_tx_intron['intron'].iloc[0] != matches['intron'].iloc[0]:
                print("multiple introns involved")
                start_intron = int(matches['intron'].iloc[0])
                end_intron = int(within_tx_intron['intron'].iloc[0])
                supp_event = f"exon {'-'.join(str(i) for i in range(start_intron+1, end_intron+1))} skipping/"
                print(f"the second event is {supp_event}")
                distance = start_end - within_tx_intron['end'].iloc[0]
                direction = "+" if distance > 0 else ""
            else:
                supp_event = ""    
                
            event = f"{supp_event}intronic cryptic {splice_site_type} @ {direction}{distance}"
            return {'transcript': tx, 'event': event}
            
        else:
            print("didn't find an intronic variant")
            within_tx_exon = self.exons[(self.exons['transcript'] == tx) &
                            (self.exons['start'] <= start_end) &
                            (self.exons['end'] >= start_end)
            ]
            if not within_tx_exon.empty:
                print("found an exonic variant")
                print(f"{start_end} is between {within_tx_exon}")
                event = f"exonic cryptic {splice_site_type} @ {direction}{distance}"
                return {'transcript': tx, 'event': event}
            else:
                print("no exonic variant either")
                return {"something has gone wrong!?"}