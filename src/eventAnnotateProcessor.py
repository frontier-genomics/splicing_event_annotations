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

            transcripts = set(unique_start_transcripts) & set(unique_end_transcripts)

            intron_match_transcripts = []

            for transcript in transcripts:
                start_intron = alternate_start_matches[alternate_start_matches['transcript'] == transcript]['intron'].iloc[0]
                end_intron = alternate_end_matches[alternate_end_matches['transcript'] == transcript]['intron'].iloc[0]
                if start_intron == end_intron:
                    intron_match_transcripts.append(transcript)

            transcripts = intron_match_transcripts

            start_intron = alternate_start_matches['intron'].iloc[0]
            end_intron = alternate_end_matches['intron'].iloc[0]
            start_tx = alternate_start_matches['transcript'].iloc[0]
            if self.coordinates['type'] == "ir":
                return {"no current support for alternate intron retention"}
            elif start_intron == end_intron:
                if len(transcripts) > 1:
                    print("all alternate start and end match transcripts match each other")
                    print(f"Number of unique start matches: {len(unique_start_transcripts)}")
                    print(f"Number of unique end matches: {len(unique_end_transcripts)}")
                    print(unique_start_transcripts)
                    #transcripts = set(unique_start_transcripts) & set(unique_end_transcripts)
                    print(transcripts)
                    transcripts_str = ', '.join(sorted(transcripts))
                    event = f"{transcripts_str}"
                    return {'alternate': 'alternate ', 'event': f" ({event})"}
                else:
                    event = f"exon {start_intron}-{end_intron+1}"
                    return {'alternate': 'alternate ', 'event': f" ({transcripts[0]} {event})"}
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
                distance = distance + 1 if self.coordinates['strand'] == "-" and splice_site_type == "acceptor" else distance
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
            location = "intronic junction"
            alternate = self._is_alternate()
            unannotated = "unannotated " if alternate['alternate'] == "" else ""
            event = f"{alternate['alternate']}{unannotated}{location}{strand_in}{alternate['event']}"
            return {'transcript': tx, 'event': event}
        
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
            








































    def reference_match(self, start_end, base=1):
        """
        Finds a reference match for the given start_end position.

        Args:
            start_end (str): The position to find a reference match for. Should be either 'start' or 'end'.
            base (int, optional): The base value to subtract from the start position if start_end is 'start'. Defaults to 1.

        Returns:
            pandas.DataFrame: A DataFrame containing the matching rows from the annotation data.

        Raises:
            ValueError: If an invalid start_end choice is provided.

        """
        print(f"finding reference match for {start_end}")
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
        return match
    


    def _produce_annotation(self, annotation, start_matches, end_matches):
        """
        Produces an annotation string based on the given annotation dictionary.

        Args:
            annotation (dict): A dictionary containing the annotation information.

        Returns:
            str: The generated annotation string.

        """
        event = annotation['event']
        cryptic = annotation['cryptic']
        supplementary_event = annotation['supplementary_event']
        location = annotation['location']
        distance = annotation['distance']
        direction = annotation['direction']

        if cryptic not in ['canonical ', 'intron ']:  # Replace 'event1', 'event2', 'event3' with the actual events
            alternate = self._is_alternate_splicing(start_matches, end_matches)['alternate']
            alternate_event = self._is_alternate_splicing(start_matches, end_matches)['alternate_event']
            if alternate == "alternate ":
                cryptic = ""
        else:
            alternate = ""
            alternate_event = ""

        print(cryptic)
        print(supplementary_event)
        print(location)
        print(event)
        print(direction)
        print(distance)
        print(alternate)
        print(alternate_event)

        return f"{alternate}{supplementary_event}{cryptic}{location}{event}{direction}{distance}{alternate_event}"



    def fetch_transcript_annotations(self, start_matches_all_tx, end_matches_all_tx):
        transcript = self.coordinates['transcript']
        print(transcript)

        start_matches = start_matches_all_tx[start_matches_all_tx['transcript'] == transcript]
        end_matches = end_matches_all_tx[end_matches_all_tx['transcript'] == transcript]

        print(f"matching starts are \n{start_matches}")
        print(f"matching ends are \n{end_matches}")

        if end_matches.empty and start_matches.empty:
            print("no start or end matches for transcript")
            annotation_dict = self._annotate_unannotated()
            return self._produce_annotation(annotation_dict, start_matches_all_tx, end_matches_all_tx)
        
        elif start_matches.empty or end_matches.empty:
            if end_matches.empty:
                print("start matches")
                annotation_dict = self._annotate_cryptic('end', self.coordinates['end'], start_matches)
                return self._produce_annotation(annotation_dict, start_matches_all_tx, end_matches_all_tx)
            
            elif start_matches.empty:
                print("end matches")
                annotation_dict = self._annotate_cryptic('start', self.coordinates['start'], end_matches)
                return self._produce_annotation(annotation_dict, start_matches_all_tx, end_matches_all_tx)

        else:
            if self.coordinates['type'] == "ir":
                print("splicing type is intron retention")
                annotation_dict = self._annotate_intron_retention(start_matches['intron'].unique())
                return self._produce_annotation(annotation_dict, start_matches_all_tx, end_matches_all_tx)

            elif start_matches['intron'].unique() == end_matches['intron'].unique():
                print("start and end introns are the same")
                annotation_dict = self._annotate_canonical(start_matches['intron'].unique())
                return self._produce_annotation(annotation_dict, start_matches_all_tx, end_matches_all_tx)
            
            elif start_matches['intron'].unique() != end_matches['intron'].unique():
                print("start and end introns are different")
                annotation_dict = self._annotate_exon_skipping(start_matches['intron'].unique(), end_matches['intron'].unique())
                return self._produce_annotation(annotation_dict, start_matches_all_tx, end_matches_all_tx)

            return({'transcript': "unknown", 'event': "unknown event"})

    def _is_alternate_splicing(self, start_matches, end_matches):
        transcript = self.coordinates['transcript']

        start_matches = start_matches[start_matches['transcript'] != transcript]
        end_matches = end_matches[end_matches['transcript'] != transcript]

        alternate_start_matches = start_matches.sort_values(by='transcript')
        alternate_end_matches = end_matches.sort_values(by='transcript')

        print(f"matching alternate starts are \n{alternate_start_matches}")
        print(f"matching alternate ends are \n{alternate_end_matches}")

        if alternate_start_matches.empty or alternate_end_matches.empty:
            return {'alternate': '', 'alternate_event': ""}
        else:
            unique_start_transcripts = alternate_start_matches['transcript'].unique()
            unique_end_transcripts = alternate_end_matches['transcript'].unique()

            transcripts = set(unique_start_transcripts) & set(unique_end_transcripts)

            intron_match_transcripts = []

            for transcript in transcripts:
                start_intron = alternate_start_matches[alternate_start_matches['transcript'] == transcript]['intron'].iloc[0]
                end_intron = alternate_end_matches[alternate_end_matches['transcript'] == transcript]['intron'].iloc[0]
                if start_intron == end_intron:
                    intron_match_transcripts.append(transcript)

            transcripts = intron_match_transcripts

            start_intron = alternate_start_matches['intron'].iloc[0]
            end_intron = alternate_end_matches['intron'].iloc[0]

            if self.coordinates['type'] == "ir":
                return {"no current support for alternate intron retention"}
            elif start_intron == end_intron:
                if len(transcripts) > 1:
                    print("all alternate start and end match transcripts match each other")
                    print(f"Number of unique start matches: {len(unique_start_transcripts)}")
                    print(f"Number of unique end matches: {len(unique_end_transcripts)}")
                    print(unique_start_transcripts)
                    #transcripts = set(unique_start_transcripts) & set(unique_end_transcripts)
                    print(transcripts)
                    transcripts_str = ', '.join(sorted(transcripts))
                    return {'alternate': 'alternate ', 'cryptic': '', 'alternate_event': f" ({transcripts_str})"}
                else:
                    event = f"exon {start_intron}-{end_intron+1}"
                    return {'alternate': 'alternate ', 'cryptic': '', 'alternate_event': f" ({transcripts[0]} {event})"}
            else:
                return {'alternate': '', 'alternate_event': ""}

    def _annotate_cryptic(self, start_end, position, matches):
        #Set cryptic to cryptic
        cryptic = "cryptic "
        print("identified a cryptic")

        #Set event to acceptor/donor and calculate distance, based on strand
        strand = self.coordinates['strand']
        start = self.coordinates['start']
        end = self.coordinates['end']

        #Set location based on within_tx_intron and within_tx_exon
        transcript = self.coordinates['transcript']

        within_tx_intron = self.annotation[(self.annotation['transcript'] == transcript) &
                        (self.annotation['start'] <= position) &
                        (self.annotation['end'] >= position)
        ]
        if start_end == "start":
            position = matches['start'].unique()[0] 
        elif start_end == "end":
            position = matches['end'].unique()[0] 

        if not within_tx_intron.empty:
            print("located within an intron")

            intron_1 = within_tx_intron['intron'].unique()
            intron_2 = matches['intron'].unique()

            if intron_1 != intron_2:
                print("but event spans multiple introns")
                supp_event = self._annotate_supplementary(intron_1, intron_2)
                position = within_tx_intron['start'].unique()[0] if start_end == "start" else within_tx_intron['end'].unique()[0]
            else:
                supp_event = ""    

            if start_end == "start":
                if strand == "+":
                    event = "donor"
                    distance = start - position + 1
                elif strand == "-":
                    event = "acceptor"
                    distance = position - start
            elif start_end == "end":
                if strand == "+":
                    event = "acceptor"
                    distance = end - position - 1
                elif strand == "-":
                    event = "donor"
                    distance = position - end
            print(event)

            #Set direction        
            direction = " @ +" if distance > 0 else " @ "
            
            location = "intronic "

            return {'event': event,
                    'cryptic': cryptic,
                    'supplementary_event': supp_event,
                    'location': location,
                    'distance': distance,
                    'direction': direction,
                    'alternate': ""}
            
        else:
            within_tx_exon = self.exons[(self.exons['transcript'] == transcript) &
                            (self.exons['start'] <= position) &
                            (self.exons['end'] >= position)
            ]

            if not within_tx_exon.empty:
                print("located within an exon")

                exon = within_tx_exon['exon'].unique()
                intron = matches['intron'].unique()

                if exon != intron and exon != intron+1 :
                    print("but event spans multiple introns")
                    supp_event = self._annotate_supplementary(exon, intron)
                    position = within_tx_exon['start'].unique()[0] if start_end == "start" else within_tx_exon['end'].unique()[0]
                else:
                    supp_event = ""   

                if start_end == "start":
                    if strand == "+":
                        event = "donor"
                        distance = start - position
                    elif strand == "-":
                        event = "acceptor"
                        distance = position - start
                elif start_end == "end":
                    if strand == "+":
                        event = "acceptor"
                        distance = end - position
                    elif strand == "-":
                        event = "donor"
                        distance = position - end
                print(event)

                #Set direction        
                direction = " @ +" if distance > 0 else " @ "
                
                location = "exonic "

                return {'event': event,
                        'cryptic': cryptic,
                        'supplementary_event': supp_event,
                        'location': location,
                        'distance': distance,
                        'direction': direction,
                        'alternate': ""}
        
            else:
                print("located outside of transcript boundaries")
                location = "intergenic "

                return {'event': "",
                        'cryptic': "",
                        'supplementary_event': "",
                        'location': location,
                        'distance': "",
                        'direction': "",
                        'alternate': ""}
            
    def _annotate_unannotated(self):
        tx = self.coordinates['transcript']

        within_tx_intron_start = self.annotation[
            (self.annotation['chrom'] == self.coordinates['chrom']) &
            (self.annotation['start'] <= self.coordinates['start']) &
            (self.annotation['transcript'] == tx) &
            (self.annotation['end'] >= self.coordinates['start'])
        ]
        
        within_tx_intron_end = self.annotation[
            (self.annotation['chrom'] == self.coordinates['chrom']) &
            (self.annotation['start'] <= self.coordinates['end']) &
            (self.annotation['transcript'] == tx) &
            (self.annotation['end'] >= self.coordinates['end'])
        ]
        
        if not within_tx_intron_start.empty and not within_tx_intron_end.empty:
            print("both start and end are within introns")
            strand_check = within_tx_intron_start['strand'].iloc[0] == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""
            return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "intronic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}
        
        elif not within_tx_intron_start.empty and within_tx_intron_end.empty:
            print("only start is within intron")
            strand_check = within_tx_intron_start['strand'].iloc[0] == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""
            within_tx_exon = self.exons[
                (self.exons['transcript'] == tx) &
                (self.exons['start'] <= self.coordinates['end']) &
                (self.exons['end'] >= self.coordinates['end'])
            ]
            print(within_tx_exon)
            if not within_tx_exon.empty:
                return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "intronic/exonic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}
            else:
                return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "intergenic/intronic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}
        
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
                return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "intronic/exonic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}
            else:
                return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "intergenic/intronic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}
            
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
                strand_check = within_tx_exon_start['strand'].iloc[0] == self.coordinates['strand']
                strand_in = " (opposite strand)" if not strand_check else ""

                return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "exonic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}
        
            elif not within_tx_exon_start.empty and within_tx_exon_end.empty:
                print("only start is within exon")

                strand_check = within_tx_exon_start['strand'].iloc[0] == self.coordinates['strand']
                strand_in = " (opposite strand)" if not strand_check else ""

                return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "intergenic/exonic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}

            # might be overkill to check for both start and end separately
            elif within_tx_exon_start.empty and not within_tx_exon_end.empty:
                print("only end is within exon")

                strand_check = within_tx_exon_end['strand'].iloc[0] == self.coordinates['strand']
                strand_in = " (opposite strand)" if not strand_check else ""    

                return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "intergenic/exonic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}
            
            else:
                return {'event': "junction",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'location': "unknown ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}

    def _annotate_intron_retention(self, start_intron):
        intron = start_intron[0]
        cryptic = "intron "
        event = f"{intron} retention"
        
        return {'event': event,
                'cryptic': cryptic,
                'supplementary_event': "",
                'location': "",
                'distance': "",
                'direction': "",
                'alternate': ""}

    def _annotate_canonical(self, start_intron):
        cryptic = "canonical "
        intron = start_intron[0]
        event = f"exon {intron}-{intron+1} splicing"

        return {'event': event,
                'cryptic': cryptic,
                'supplementary_event': "",
                'location': "",
                'distance': "",
                'direction': "",
                'alternate': ""}

    def _annotate_exon_skipping(self, start_intron, end_intron):
        introns = [start_intron[0], end_intron[0]]
        event = f"exon {'-'.join(str(i) for i in range(min(introns)+1, max(introns)+1))} skipping"

        

        return {'event': event,
                'cryptic': "",
                'supplementary_event': "",
                'location': "",
                'distance': "",
                'direction': "",
                'alternate': ""}
    
    def _annotate_supplementary(self, start_intron, end_intron):
        introns = [start_intron[0], end_intron[0]]
        event = f"exon {'-'.join(str(i) for i in range(min(introns)+1, max(introns)+1))} skipping/"

        return event