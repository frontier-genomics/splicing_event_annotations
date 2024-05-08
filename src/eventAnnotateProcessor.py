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
        if self.annotation_choice == "refseq":
            self.annotation = pd.read_csv("resources/annotations/refseq_curated_introns_sorted.tsv", sep='\t')
            self.exons = pd.read_csv("resources/annotations/refseq_curated_exons_sorted.tsv", sep='\t')
        elif self.annotation_choice == "ensembl":
            raise ValueError("Ensembl annotations are currently not supported. Please try again with 'refseq'.")
        else:
            raise ValueError("Invalid annotation choice. Please select either 'refseq' or 'ensembl' (currently not supported).")
        
        self.annotation.columns = ['chrom', 'start', 'end', 'transcript', 'intron', 'strand']

    def get_mane_transcript(self, base=1):
        print(f"finding mane transcript match for event")
        transcripts = pd.read_csv("resources/annotations/refseq_mane_gene_tx_names.tsv", sep='\t')

        chrom = self.coordinates['chrom']
        start = self.coordinates['start']
        end = self.coordinates['start']
        strand = self.coordinates['strand']

        matching_rows_start = transcripts[
            (transcripts['chrom'] == chrom) &
            (transcripts['start'] <= start) &
            (transcripts['end'] >= start)
        ]

        matching_rows_end = transcripts[
            (transcripts['chrom'] == chrom) &
            (transcripts['start'] <= end) &
            (transcripts['end'] >= end)
        ]

        matching_rows = pd.concat([matching_rows_start, matching_rows_end]).drop_duplicates()

        if matching_rows.empty:
            print("no MANE transcript match found")
            return {"transcript": "unknown",
                    "gene": "unknown",
                    "warning": "no MANE transcript match found"}
        
        matching_rows['start'] = matching_rows['start'] + base
        match = pd.DataFrame(matching_rows)

        print(match)
        print(len(match['transcript'].unique()))

        if len(match['transcript'].unique()) == 1:
            transcript = match['transcript'].iloc[0]
            gene = match['gene'].iloc[0]

            if match['strand'].unique()[0] == strand:

                warning = "none"
            else:
                print(match['strand'].unique()[0])
                print(strand)
                warning = "MANE transcript found on opposite strand only."

        elif len(match['transcript'].unique()) > 1:
            if sum(match['strand'] == strand) == 1:
                transcript = match['transcript'].iloc[0]
                gene = match['gene'].iloc[0]
                warning = "Overlapping MANE transcript on opposite strand found."
            elif sum(match['strand'] == strand) != 1:
                transcript = "unknown"
                gene = "unknown"
                warning = "Multiple overlapping MANE transcripts found. Unable to assign."
            elif sum(match['strand'] == strand) == 0:
                transcript = "unknown"
                gene = "unknown"
                warning = "Multiple overlapping MANE transcripts found on opposite strand only. Unable to assign."

        return {"transcript": transcript,
                "gene": gene,
                "warning": warning}

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
        introns = annotation['introns']
        event_type = annotation['event_type']

        if cryptic not in ['canonical ', 'intron ']:
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

        return {'event': f"{alternate}{supplementary_event}{cryptic}{location}{event}{direction}{distance}{alternate_event}",
                'event_type': event_type,
                'introns': introns}

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
            print(f"Matching Transcripts {transcripts}")

            intron_match_transcripts = []

            for transcript in transcripts:
                start_intron = alternate_start_matches[alternate_start_matches['transcript'] == transcript]['intron'].iloc[0]
                end_intron = alternate_end_matches[alternate_end_matches['transcript'] == transcript]['intron'].iloc[0]
                if start_intron == end_intron:
                    intron_match_transcripts.append(transcript)

            transcripts = intron_match_transcripts
            print(f"Matching Transcripts After Filtering {transcripts}")

            start_intron = alternate_start_matches['intron'].iloc[0]
            end_intron = alternate_end_matches['intron'].iloc[0]

            if self.coordinates['type'] == "ir":
                return {"no current support for alternate intron retention"}
            else:
                if len(transcripts) > 1:
                    print("all alternate start and end match transcripts match each other")
                    print(f"Number of unique start matches: {len(unique_start_transcripts)}")
                    print(f"Number of unique end matches: {len(unique_end_transcripts)}")
                    print(unique_start_transcripts)
                    print(transcripts)
                    transcripts_str = ', '.join(sorted(transcripts))
                    return {'alternate': 'alternate ', 'cryptic': '', 'alternate_event': f" ({transcripts_str})"}
                elif len(transcripts) == 1:
                    event = f"exon {start_intron}-{end_intron+1}"
                    return {'alternate': 'alternate ', 'cryptic': '', 'alternate_event': f" ({transcripts[0]} {event})"}
                else:
                    return {'alternate': '', 'alternate_event': ""}

    def _annotate_cryptic(self, start_end, var_position, matches):
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
                        (self.annotation['start'] <= var_position) &
                        (self.annotation['end'] >= var_position)
        ]
        if start_end == "start":
            position = matches['start'].unique()[0]
        elif start_end == "end":
            position = matches['end'].unique()[0] 

        if not within_tx_intron.empty:
            print("located within an intron")

            position = within_tx_intron['start'].unique()[0] if start_end == "start" else within_tx_intron['end'].unique()[0]

            introns = [int(within_tx_intron['intron'].unique()), int(matches['intron'].unique())]

            location = f"intron {int(introns[0])} "

            if introns[0] != introns[1]:
                print("but event spans multiple introns")
                supp = self._annotate_supplementary(introns)
                supp_event = supp['event']
                supp_event_type = supp["supp_event_type"]
                introns[0] = supp['introns']
            else:
                supp_event = ""
                supp_event_type = ""    

            if start_end == "start":
                if strand == "+":
                    event = "donor"
                    distance = start - position
                    print(f"The calculation is {start} - {position} = {distance}")
                elif strand == "-":
                    event = "acceptor"
                    distance = position - start
                    print(f"The calculation is {position} - {start} = {distance}")
            elif start_end == "end":
                if strand == "+":
                    event = "acceptor"
                    distance = end - position - 1
                    print(f"The calculation is {end} - {position} - 1 = {distance}")
                elif strand == "-":
                    event = "donor"
                    distance = position - end + 1
                    print(f"The calculation is {position} - {end} + 1 = {distance}")
            print(event)

            #Set direction        
            direction = " @ +" if distance > 0 else " @ "

            return {'event': event,
                    'cryptic': cryptic,
                    'supplementary_event': supp_event,
                    'location': location,
                    'distance': distance,
                    'direction': direction,
                    'alternate': "",
                    'event_type': f"{supp_event_type}{cryptic}{event}",
                    'introns': introns[0]}
            
        else:
            within_tx_exon = self.exons[(self.exons['transcript'] == transcript) &
                            (self.exons['start'] <= var_position) &
                            (self.exons['end'] >= var_position)
            ]

            if not within_tx_exon.empty:
                print("located within an exon")

                position = within_tx_exon['start'].unique()[0] if start_end == "end" else within_tx_exon['end'].unique()[0]

                introns = [int(within_tx_exon['exon'].unique()), int(matches['intron'].unique())]

                if start_end == "start":
                    if strand == "+":
                        event = "donor"
                        distance = start - position - 1
                        print(f"The calculation is {start} - {position} - 1 = {distance}")
                    elif strand == "-":
                        event = "acceptor"
                        distance = position - start + 1
                        print(f"The calculation is {position} - {start} + 1= {distance}")
                elif start_end == "end":
                    if strand == "+":
                        event = "acceptor"
                        distance = end - position
                        print(f"The calculation is {end} - {position} = {distance}")
                    elif strand == "-":
                        event = "donor"
                        distance = position - end
                        print(f"The calculation is {position} - {end} = {distance}")

                print(event)

                location = f"exon {introns[0]} "

                if introns[0] != introns[1] and introns[0] != introns[1]+1 :
                    print("but event spans multiple introns")
                    if event == 'acceptor':
                        introns[0] = introns[0] - 1
                    supp = self._annotate_supplementary(introns)
                    supp_event = supp['event']
                    supp_event_type = supp["supp_event_type"]
                    introns[0] = supp['introns']

                else:
                    supp_event = ""   
                    supp_event_type = ""
                    if event == 'acceptor':
                        introns[0] = introns[0] - 1

                #Set direction        
                direction = " @ +" if distance > 0 else " @ "

                return {'event': event,
                        'cryptic': cryptic,
                        'supplementary_event': supp_event,
                        'location': location,
                        'distance': distance,
                        'direction': direction,
                        'alternate': "",
                        'event_type': f"{supp_event_type}{cryptic}{event}",
                        'introns': introns[0]}
        
            else:
                print("located outside of transcript boundaries")
                location = "intergenic "

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

                direction = " @ +" if distance > 0 else " @ "

                return {'event': event,
                        'cryptic': "",
                        'supplementary_event': "",
                        'location': location,
                        'distance': distance,
                        'direction': direction,
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
                    'location': "intergenic ",
                    'distance': "",
                    'direction': "",
                    'alternate': ""}

    def _annotate_intron_retention(self, start_intron):
        intron = start_intron[0]
        cryptic = "intron "
        event = f"{intron} retention"
        event_type = "intron retention"

        return {'event': event,
                'cryptic': cryptic,
                'supplementary_event': "",
                'location': "",
                'distance': "",
                'direction': "",
                'alternate': "",
                'event_type': event_type,
                'introns': intron}

    def _annotate_canonical(self, start_intron):
        cryptic = "canonical "
        intron = start_intron[0]
        event = f"exon {intron}-{intron+1} splicing"
        event_type = "canonical"

        return {'event': event,
                'cryptic': cryptic,
                'supplementary_event': "",
                'location': "",
                'distance': "",
                'direction': "",
                'alternate': "",
                'event_type': event_type,
                'introns': intron}

    def _annotate_exon_skipping(self, start_intron, end_intron):
        introns = [start_intron[0], end_intron[0]]
        event = f"exon {'-'.join(str(i) for i in range(min(introns)+1, max(introns)+1))} skipping"
        event_type = "exon skipping"
        introns = f"{', '.join(str(i) for i in range(min(introns), max(introns)+1))}"

        

        return {'event': event,
                'cryptic': "",
                'supplementary_event': "",
                'location': "",
                'distance': "",
                'direction': "",
                'alternate': "",
                'event_type': event_type,
                'introns': introns}
    
    def _annotate_supplementary(self, introns):
        event = f"exon {'-'.join(str(i) for i in range(min(introns)+1, max(introns)+1))} skipping/"
        introns = f"{', '.join(str(i) for i in range(min(introns), max(introns)+1))}"

        return {'event': event,
                'supp_event_type': 'exon skipping, ',
                'introns': introns}