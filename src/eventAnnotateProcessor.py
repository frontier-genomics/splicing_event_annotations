import pandas as pd
import numpy as np

class EventAnnotate:

    def __init__(self, chrom, start, end, strand, transcript, type, annotation_choice):
        self.coordinates = {
            'chrom': str(chrom),
            'start': np.int64(start),
            'end': np.int64(end),
            'strand': str(strand),
            'transcript': str(transcript),
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

    def matcher(self, start_end, alternative = False, base = 1):
        if start_end not in ['start', 'end']:
            raise ValueError("Invalid start_end choice. Please select either 'start' or 'end'.")
        
        if start_end == "start":
            query = self.coordinates[start_end] - base
        elif start_end == "end":
            query = self.coordinates[start_end]
        
        chrom = self.coordinates['chrom']
        transcript = self.coordinates['transcript'] 
        
        if alternative:
            matching_rows = self.annotation[
                (self.annotation[start_end] == query) &
                (self.annotation['chrom'] == chrom) &
                (self.annotation['transcript'] != transcript)
            ]
        else:
            matching_rows = self.annotation[
                (self.annotation[start_end] == query) &
                (self.annotation['chrom'] == chrom) &
                (self.annotation['transcript'] == transcript)
            ]
        
        matching_rows['start'] = matching_rows['start'] + base
        match = pd.DataFrame(matching_rows)
        return(match)
        
    def create_annotations(self, start_matches, end_matches):
        if end_matches.empty and start_matches.empty:
            #check if either start or end matches are within the annotation for the transcript
            return self._annotate_unannotated()
        
        elif start_matches.empty or end_matches.empty:
            if end_matches.empty:
                print("start matches")
                end_matches = ""
                
                start_intron = start_matches['intron'].iloc[0]
                start_tx = start_matches['transcript'].iloc[0]
                
                match_coords = start_matches['end'].iloc[0]
                self_coords = self.coordinates['end']
                
                strand = self.coordinates['strand']
                splice_site_type = "acceptor" if strand == "+" else "donor"
                distance = self_coords - match_coords if strand == "+" else match_coords - self_coords
                
                event_name = self._name_event(start_tx, distance, self.coordinates['end'], splice_site_type, start_matches)

                return(event_name)
            
            elif start_matches.empty:
                print("end matches")
                start_matches = ""
                
                end_intron = end_matches['intron'].iloc[0]
                end_tx = end_matches['transcript'].iloc[0]
                
                match_coords = end_matches['start'].iloc[0]
                self_coords = self.coordinates['start']
                
                strand = self.coordinates['strand']
                splice_site_type = "donor" if strand == "+" else "acceptor"
                distance = self_coords - match_coords if strand == "+" else match_coords - self_coords
    
                event_name = self._name_event(end_tx, distance, self.coordinates['start'], splice_site_type, end_matches)

                return(event_name)
        
        else:
            start_intron = start_matches['intron'].iloc[0]
            end_intron = end_matches['intron'].iloc[0]
            start_tx = start_matches['transcript'].iloc[0]
            end_tx = start_matches['transcript'].iloc[0]
            introns = [start_intron, end_intron]

            if start_tx == end_tx:
                if self.coordinates['type'] == "ir":
                    event = f"intron {start_intron} retention"
                    return {'transcript': start_tx, 'event': event}
                elif start_intron == end_intron:
                    event = f"canonical exon {start_intron}-{end_intron+1} splicing"
                    return {'transcript': start_tx, 'event': event}
                elif start_intron != end_intron:
                    event = f"exon {'-'.join(str(i) for i in range(min(introns)+1, max(introns)+1))} skipping"
                    alternate = self._is_alternate()
                    event = f"{alternate['alternate']}{event}{alternate['event']}"
                    return {'transcript': start_tx, 'event': event}

            elif start_tx != end_tx:
                    return {'transcript': "transcript mismatch", 'event': "unknown event"}

    def _is_alternate(self):
        alternate_start_matches = self.matcher('start', alternative = True)
        alternate_end_matches = self.matcher('end', alternative = True)
        alternate_start_matches = alternate_start_matches.sort_values(by='transcript')
        alternate_end_matches = alternate_end_matches.sort_values(by='transcript')
        print(alternate_start_matches)
        print(alternate_end_matches)

        if alternate_start_matches.empty or alternate_end_matches.empty:
            return {'alternate': "", 'event': ""}
        else:
            unique_start_transcripts = alternate_start_matches['transcript'].unique()
            unique_end_transcripts = alternate_end_matches['transcript'].unique()
            start_intron = alternate_start_matches['intron'].iloc[0]
            end_intron = alternate_end_matches['intron'].iloc[0]
            start_tx = alternate_start_matches['transcript'].iloc[0]
            if self.coordinates['type'] == "ir":
                event = f"intron {start_intron} retention"
                return {'alternate': 'alternate ', 'event': f" ({start_tx} {event})"}
            elif start_intron == end_intron:
                if len(unique_start_transcripts) > 1 and len(unique_end_transcripts) > 1:
                    print("all alternate start and end match transcripts match each other")
                    print(f"Number of unique start matches: {len(unique_start_transcripts)}")
                    print(f"Number of unique end matches: {len(unique_end_transcripts)}")
                    print(unique_start_transcripts)
                    transcripts = set(unique_start_transcripts) & set(unique_end_transcripts)
                    print(transcripts)
                    transcripts_str = ', '.join(sorted(transcripts))
                    event = f"{transcripts_str}"
                    return {'alternate': 'alternate ', 'event': f" ({event})"}
                else:
                    event = f"exon {start_intron}-{end_intron+1}"
                    return {'alternate': 'alternate ', 'event': f" ({start_tx} {event})"}
            else:
                return {'alternate': "", 'event': ""}

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
                introns = [start_intron, end_intron]
                supp_event = f"exon {'-'.join(str(i) for i in range(min(introns)+1, max(introns)+1))} skipping/"
                print(f"the second event is {supp_event}")
                distance = start_end - within_tx_intron['end'].iloc[0]
                direction = "+" if distance > 0 else ""
            else:
                supp_event = ""    
                
            distance = distance - 1 if self.coordinates['strand'] == "+" and splice_site_type == "acceptor" else distance
            distance = distance + 1 if self.coordinates['strand'] == "-" and splice_site_type == "acceptor" else distance
            distance = distance + 1 if self.coordinates['strand'] == "+" and splice_site_type == "donor" else distance
            
            location = "intronic "
            alternate = self._is_alternate()
            print(alternate)
            cryptic = "cryptic " if alternate['alternate'] == "" else ""
            event = f"{alternate['alternate']}{supp_event}{location}{cryptic}{splice_site_type} @ {direction}{distance}{alternate['event']}"
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

                if within_tx_exon['exon'].iloc[0] != matches['intron'].iloc[0] and within_tx_exon['exon'].iloc[0] != matches['intron'].iloc[0]+1 :
                    print("multiple exons/introns involved")
                    print(f"match with intron {matches['intron'].iloc[0]}")
                    print(f"within exon {within_tx_exon['exon'].iloc[0]}")
                    print(f"exon != intron {within_tx_exon['exon'].iloc[0] != matches['intron'].iloc[0]}")
                    print(f"exon != intron+1 {within_tx_exon['exon'].iloc[0] != matches['intron'].iloc[0]+1}")
                    start_intron = int(matches['intron'].iloc[0])
                    end_intron = int(within_tx_exon['exon'].iloc[0])
                    introns = [start_intron, end_intron]
                    supp_event = f"exon {'-'.join(str(i) for i in range(min(introns)+1, max(introns)+1))} skipping/"
                    print(f"the second event is {supp_event}")
                    distance = start_end - within_tx_exon['end'].iloc[0] - 1
                    direction = "+" if distance > 0 else ""
                else:
                    supp_event = ""   

                distance = distance + 1 if self.coordinates['strand'] == "-" and splice_site_type == "donor" else distance
                
                event = f"{supp_event}exonic cryptic {splice_site_type} @ {direction}{distance}"
                return {'transcript': tx, 'event': event}
            else:
                print("no exonic variant either")
                location = "intergenic "
                alternate = self._is_alternate()
                cryptic = "cryptic " if alternate['alternate'] == "" else ""
                event = f"{alternate['alternate']}{location}{cryptic}{splice_site_type} @ {direction}{distance}{alternate['event']}"
                return {'transcript': tx, 'event': event}
                # Find alternative splicing
                # Find unannotated splicing
            
    def _annotate_unannotated(self):
        within_tx_intron_start = self.annotation[
            (self.annotation['chrom'] == self.coordinates['chrom']) &
            (self.annotation['start'] <= self.coordinates['start']) &
            (self.annotation['transcript'] == self.coordinates['transcript']) &
            (self.annotation['end'] >= self.coordinates['start'])
        ]
        
        within_tx_intron_end = self.annotation[
            (self.annotation['chrom'] == self.coordinates['chrom']) &
            (self.annotation['start'] <= self.coordinates['end']) &
            (self.annotation['transcript'] == self.coordinates['transcript']) &
            (self.annotation['end'] >= self.coordinates['end'])
        ]
        
        if not within_tx_intron_start.empty and not within_tx_intron_end.empty:
            print("both start and end are within introns")
            tx = within_tx_intron_start['transcript'].iloc[0]
            strand_check = within_tx_intron_start['strand'].iloc[0] == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""
            return {'transcript': tx, 'event': f"unannotated intronic junction{strand_in}"}
        
        elif not within_tx_intron_start.empty and within_tx_intron_end.empty:
            print("only start is within intron")
            tx = within_tx_intron_start['transcript'].iloc[0]
            strand_check = within_tx_intron_start['strand'].iloc[0] == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""
            print(tx)
            within_tx_exon = self.exons[
                (self.exons['transcript'] == tx) &
                (self.exons['start'] <= self.coordinates['end']) &
                (self.exons['end'] >= self.coordinates['end'])
            ]
            print(within_tx_exon)
            if not within_tx_exon.empty:
                return {'transcript': tx, 'event': f"unannotated intronic/exonic junction{strand_in}"}
            else:
                return {'transcript': tx, 'event': f"unknown intergenic/intronic junction{strand_in}"}
        
        elif within_tx_intron_start.empty and not within_tx_intron_end.empty:
            print("only end is within intron")
            tx = within_tx_intron_end['transcript'].iloc[0]
            strand_check = within_tx_intron_end['strand'].iloc[0] == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""    
            
                        
            within_tx_exon = self.exons[
                (self.exons['transcript'] == tx) &
                (self.exons['start'] <= self.coordinates['start']) &
                (self.exons['end'] >= self.coordinates['start'])
            ]
            
            if not within_tx_exon.empty:
                return {'transcript': tx, 'event': f"unannotated intronic/exonic junction{strand_in}"}
            else:
                return {'transcript': tx, 'event': f"unknown intergenic/intronic junction{strand_in}"}
            
        else:
            within_tx_exon_start = self.exons[
                (self.exons['chr'] == self.coordinates['chrom']) &
                (self.exons['start'] <= self.coordinates['start']) &
                (self.exons['transcript'] == self.coordinates['transcript']) &
                (self.exons['end'] >= self.coordinates['start'])
            ]

            within_tx_exon_end = self.exons[
                (self.exons['chr'] == self.coordinates['chrom']) &
                (self.exons['start'] <= self.coordinates['end']) &
                (self.exons['transcript'] == self.coordinates['transcript']) &
                (self.exons['end'] >= self.coordinates['end'])
            ]
            
            if not within_tx_exon_start.empty and not within_tx_exon_end.empty:
                print("both start and end are within exons")
                tx = within_tx_exon_start['transcript'].iloc[0]
                strand_check = within_tx_exon_start['strand'].iloc[0] == self.coordinates['strand']
                strand_in = " (opposite strand)" if not strand_check else ""
                return {'transcript': tx, 'event': f"unannotated exonic junction{strand_in}"}
        
            elif not within_tx_exon_start.empty and within_tx_exon_end.empty:
                print("only start is within exon")
                tx = within_tx_exon_start['transcript'].iloc[0]
                strand_check = within_tx_exon_start['strand'].iloc[0] == self.coordinates['strand']
                strand_in = " (opposite strand)" if not strand_check else ""
                return {'transcript': tx, 'event': f"unknown intergenic/exonic junction{strand_in}"}

            elif within_tx_exon_start.empty and not within_tx_exon_end.empty:
                print("only end is within exon")
                tx = within_tx_exon_end['transcript'].iloc[0]
                strand_check = within_tx_exon_end['strand'].iloc[0] == self.coordinates['strand']
                strand_in = " (opposite strand)" if not strand_check else ""    
                return {'transcript': tx, 'event': f"unknown intergenic/exonic junction{strand_in}"}
            
            else:
                return {'transcript': "unknown", 'event': f"unknown junction"}