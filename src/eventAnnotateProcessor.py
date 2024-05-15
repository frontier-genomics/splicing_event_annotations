import pandas as pd
import numpy as np
import logging

logging.basicConfig(level=logging.INFO)

class EventAnnotate:

    def __init__(self, chrom, start, end, strand, transcript, type, dataset):
        self.coordinates = {
            'chrom': str(chrom),
            'start': int(start),
            'end': int(end),
            'strand': str(strand),
            'transcript': str(transcript),
            'type': str(type)
        }

        print(f"the current transcript for this event is {self.coordinates['transcript']}")

        self.read_refgene(dataset)

        if self.coordinates['transcript'] == "NA":
            self.coordinates['transcript'] = self.get_mane_transcript()['transcript']

        start = self.reference_match('start')
        end = self.reference_match('end')



        annotations = self.fetch_transcript_annotations(start, end)
        self.event = annotations['event']
        self.event_type = annotations['event_type']
        self.introns = annotations['introns']
        self.location = annotations['location']
        self.distance_from_authentic = annotations['distance_from_authentic']

    # def get_annotations(self, annotation_choice):
    #     if annotation_choice == "refseq":
    #         self.annotation = pd.read_csv("resources/annotations/refseq_curated_introns_sorted.tsv", sep='\t')
    #         self.exons = pd.read_csv("resources/annotations/refseq_curated_exons_sorted.tsv", sep='\t')
    #     elif annotation_choice == "ensembl":
    #         raise ValueError("Ensembl annotations are currently not supported. Please try again with 'refseq'.")
    #     else:
    #         raise ValueError("Invalid annotation choice. Please select either 'refseq' or 'ensembl' (currently not supported).")
        
    #     self.annotation.columns = ['chrom', 'start', 'end', 'transcript', 'intron', 'strand']

    def get_mane_transcript(self, base=1):
        print(f"finding mane transcript match for event")
        #transcripts = pd.read_csv("resources/annotations/refseq_mane_gene_tx_names.tsv", sep='\t')

        chrom = self.coordinates['chrom']
        start = self.coordinates['start']
        end = self.coordinates['start']
        strand = self.coordinates['strand']

        print(f"{chrom} {start} {end} {strand}")

        #print([entry for entry in self.refgene if entry['chrom'] == chrom and entry['start'] <= start and entry['end'] >= start and entry['mane'] == "Y"])
        #print([entry for entry in self.refgene if entry['chrom'] == chrom and entry['start'] <= end and entry['end'] >= end and entry['mane'] == "Y"])

        # matching_rows_start = transcripts[
        #     (transcripts['chrom'] == chrom) &
        #     (transcripts['start'] <= start) &
        #     (transcripts['end'] >= start)
        # ]

        # matching_rows_end = transcripts[
        #     (transcripts['chrom'] == chrom) &
        #     (transcripts['start'] <= end) &
        #     (transcripts['end'] >= end)
        # ]

        matching_rows_start = [entry for entry in self.refgene if entry['start'] <= start and entry['end'] >= start and entry['chrom'] == chrom and entry['mane'] == 'Y']
        matching_rows_end = [entry for entry in self.refgene if entry['start'] <= end and entry['end'] >= end and entry['chrom'] == chrom and entry['mane'] == 'Y']
        
        print(f"matching start: \n {matching_rows_start}")
        print(f"matching end: \n {matching_rows_end}")

        matching_rows = matching_rows_start + matching_rows_end


        matching_entries = {}
        for entry in matching_rows:
            matching_entries[entry['id']] = [entry['id'], entry['gene_name'], entry['strand']]

        matching_entries_unique = [entry for entry in matching_entries.values()]

        print(f"matching transcripts: \n {matching_entries_unique}")

        if not matching_entries_unique:
            print("no MANE transcript match found")
            return {"transcript": "unknown",
                    "gene": "unknown",
                    "warning": "no MANE transcript match found"}

        if len(matching_entries_unique) == 1:
            transcript = matching_entries_unique[0]

            if transcript[2] == strand:
                warning = "none"
            else:
                print(transcript[2] == strand)
                print(strand)
                warning = "MANE transcript found on opposite strand only."

        elif len(matching_entries_unique) > 1:
            strand_matches = []
            other_strand_matches = []
            for match in matching_entries_unique:
                if match[2] == strand:
                    strand_matches.append(match)
                elif match[2] != strand:
                    other_strand_matches.append(match)
            if len(strand_matches) == 1:
                transcript = strand_matches[0]
                warning = "MANE transcript on opposite strand also present."
            else:
                if len(strand_matches) > 1:
                    transcript = "NA"
                    gene = "NA"
                    warning = "Multiple MANE transcripts found on same strand. Unable to assign transcript. Please supply."
                elif len(other_strand_matches) > 1:
                    transcript = "NA"
                    gene = "NA"
                    warning = "Multiple MANE transcripts found on opposite strand. Unable to assign transcript. Please supply."


        print(f"MANE match: {transcript[0]}")
        return {"transcript": transcript[0],
                "gene": transcript[1],
                "warning": warning}

    def reference_match(self, start_end):
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
            query = self.coordinates[start_end]
        elif start_end == "end":
            query = self.coordinates[start_end]

        chrom = self.coordinates['chrom']

        print(query)

        matching_rows = [entry for entry in self.refgene if entry['start'] <= query and entry['end'] >= query and entry['chrom'] == chrom]

        region_information = self.get_exon_intron(matching_rows, query)

        return {'matching_rows': matching_rows,
                'region_information': region_information
                }
    
    def get_exon_intron(self, matching_rows, query, base = 0):
        region_dict = {}
        for transcript in matching_rows:
            tx_id = transcript['id']
            indices = transcript['exons']
            strand = transcript['strand']
            
            no_of_exons = len(indices)

            for i in range(len(indices)-1):
                if indices[i] <= query <= indices[i+1]:
                    index_match = i
                    index_position_start = indices[i]
                    index_position_end = indices[i+1]
                    match = False
                    break

            if index_match % 2 != 0:
                region_type = "intron"
                match = query == index_position_start+1 or query == index_position_end
            else:
                region_type = "exon"

            if strand == "-":
                region_number = int((no_of_exons - index_match) / 2).__floor__()
                print(f"\tthe region is {region_type} {region_number}")
            else:
                region_number = int((index_match) / 2 + 1).__floor__()
                print(f"\tthe region is {region_type} {region_number}")

            region_dict[tx_id] = {'transcript': tx_id,
                                  'region_type': region_type,
                                  'region_number': region_number,
                                  'match': match,
                                  'start-index': index_match,
                                  'end-index': index_match+1,
                                  'strand': strand}

        return region_dict

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
        location_raw = location
        distance = annotation['distance']
        distance_raw = distance
        direction = annotation['direction']
        introns = annotation['introns']
        event_type = annotation['event_type']
        supp_event_type = annotation['supp_event_type']

        if cryptic not in ['canonical ', 'intron ']:
            alternate_results = self._is_alternate_splicing(start_matches, end_matches)
            alternate = alternate_results['alternate']
            alternate_event = alternate_results['alternate_event']
            if alternate == "alternate ":
                cryptic = ""
            event_type = f"{alternate}{supp_event_type}{cryptic}{event_type}"

        else:
            alternate = ""
            alternate_event = ""

        if distance == "NA":
            distance_raw = "NA"
            distance = ""

        if location == "NA":
            location_raw = "NA"
            location = ""

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
                'introns': introns,
                'location': location_raw.strip(),
                'distance_from_authentic': distance_raw}

    def fetch_transcript_annotations(self, start_matches_all_tx, end_matches_all_tx):
        transcript = self.coordinates['transcript']
        print(transcript)

        start_matches_all = start_matches_all_tx['matching_rows']
        end_matches_all = end_matches_all_tx['matching_rows']

        start_matches = [entry for entry in start_matches_all_tx['matching_rows'] if entry['id'] == transcript]
        end_matches = [entry for entry in end_matches_all_tx['matching_rows'] if entry['id'] == transcript]

        start_region_all = start_matches_all_tx['region_information']
        end_region_all = end_matches_all_tx['region_information']

        if transcript not in list(start_region_all.keys()) and transcript not in list(end_region_all.keys()):
            print(transcript)
            print(list(start_matches_all_tx.keys()))
            print(list(end_matches_all_tx.keys()))
            print("no start or end matches for transcript")
            annotation_dict = self._annotate_unannotated(start_region_all, end_region_all)
            return self._produce_annotation(annotation_dict, start_region_all, end_region_all)
        
        if transcript not in list(start_region_all.keys()):
            start_region_all[transcript] = {'transcript': transcript,
                                  'region_type': 'intergenic',
                                  'region_number': -1,
                                  'match': False,
                                  'start-index': -1,
                                  'end-index': -1,
                                  'strand': "NA"}
        
        if transcript not in list(end_region_all.keys()):
            end_region_all[transcript] = {'transcript': transcript,
                                  'region_type': 'intergenic',
                                  'region_number': -1,
                                  'match': False,
                                  'start-index': -1,
                                  'end-index': -1,
                                  'strand': "NA"}

        start_region = start_region_all[transcript]
        end_region = end_region_all[transcript]

        start_match = start_region['match']
        end_match = end_region['match']

        start_region_all = start_matches_all_tx['region_information']
        end_region_all = end_matches_all_tx['region_information']

        # print(f"Junction start in {start_matches[0]['id']} in {start_region['region_type']} {start_region['region_number']}")
        # print(f"Junction end in {end_matches[0]['id']} in {end_region['region_type']} {end_region['region_number']}")

        if start_match:
            print("The start is annotated")
        else:
            print("The start is not annotated")

        if end_match:
            print("The end is annotated")
        else:
            print("The end is not annotated")

        print(f"matching starts are \n{start_matches}")
        print(f"matching ends are \n{end_matches}")

        print(f"start region is {start_region}")
        print(f"end region is {end_region}")


        if end_match == False and start_match == False:
            print("no start or end matches for transcript")
            annotation_dict = self._annotate_unannotated(start_region_all, end_region_all)
            return self._produce_annotation(annotation_dict, start_region_all, end_region_all)
        
        elif end_match == False or start_match == False:
            if end_match == False:
                print("start matches")
                annotation_dict = self._annotate_cryptic('end', self.coordinates['end'], start_matches, start_region_all, end_region_all)
                return self._produce_annotation(annotation_dict, start_region_all, end_region_all)
            
            elif start_match == False:
                print("end matches")
                annotation_dict = self._annotate_cryptic('start', self.coordinates['start'], end_matches, start_region_all, end_region_all)
                return self._produce_annotation(annotation_dict, start_region_all, end_region_all)

        else:
            if self.coordinates['type'] == "ir":
                print("splicing type is intron retention")
                annotation_dict = self._annotate_intron_retention(start_region)
                return self._produce_annotation(annotation_dict, start_region_all, end_region_all)

            elif start_region['region_number'] == end_region['region_number']:
                print("start and end introns are the same")
                annotation_dict = self._annotate_canonical(start_region)
                return self._produce_annotation(annotation_dict, start_region_all, end_region_all)
            
            elif start_region['region_number'] != end_region['region_number']:
                print("start and end introns are different")
                annotation_dict = self._annotate_exon_skipping(start_region, end_region)
                return self._produce_annotation(annotation_dict, start_region_all, end_region_all)

            return({'transcript': "unknown", 'event': "unknown event"})

    def _is_alternate_splicing(self, start_matches, end_matches):
        transcript = self.coordinates['transcript']

        print(start_matches)
        print(end_matches)

        # start_matches = [entry for entry in start_matches if entry['id'] != transcript]
        # end_matches = [entry for entry in end_matches if entry['id'] != transcript]

        # return {'alternate': '', 'alternate_event': ""}

        start_match = {k: v for k, v in start_matches.items() if v['transcript'] != transcript and v['match'] == True}
        end_match = {k: v for k, v in end_matches.items() if v['transcript'] != transcript and v['match'] == True}

        print(transcript)
        print(start_match)
        print(end_match)
        print(list(start_match.keys()))
        print(list(end_match.keys()))

        print("why you no print")

        filtered_data = {k: v for k, v in start_matches.items() if v['transcript'] != transcript}

        print(filtered_data)

        alternate_start_matches = start_match
        alternate_end_matches = end_match

        print(f"matching alternate starts are \n{start_matches}")
        print(f"matching alternate ends are \n{end_matches}")

        if not alternate_start_matches or not alternate_end_matches:
            return {'alternate': '', 'alternate_event': ""}
        else:
            unique_start_transcripts = list(start_match.keys())
            unique_end_transcripts = list(end_match.keys())

            transcripts = set(unique_start_transcripts) & set(unique_end_transcripts)
            print(f"Matching Transcripts {transcripts}")

            intron_match_transcripts = []

            for transcript in transcripts:
                start_intron = alternate_start_matches[transcript]['region_number']
                end_intron = alternate_end_matches[transcript]['region_number']
                if start_intron == end_intron:
                    intron_match_transcripts.append(transcript)

            transcripts = intron_match_transcripts
            print(f"Matching Transcripts After Filtering {transcripts}")

            # start_intron = alternate_start_matches['intron'].iloc[0]
            # end_intron = alternate_end_matches['intron'].iloc[0]

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
                    start_intron = alternate_start_matches[transcripts[0]]['region_number']
                    end_intron = alternate_end_matches[transcripts[0]]['region_number']
                    event = f"exon {start_intron}-{end_intron+1}"
                    return {'alternate': 'alternate ', 'cryptic': '', 'alternate_event': f" ({transcripts[0]} {event})"}
                else:
                    return {'alternate': '', 'alternate_event': ""}

    def _annotate_cryptic(self, start_end, var_position, all_matches, start_region, end_region):
        #Set cryptic to cryptic
        cryptic = "cryptic "
        print("identified a cryptic")

        #Set event to acceptor/donor and calculate distance, based on strand
        strand = self.coordinates['strand']
        query_start = self.coordinates['start']
        query_end = self.coordinates['end']

        #Set location based on within_tx_intron and within_tx_exon
        transcript = self.coordinates['transcript']

        matches = [entry for entry in all_matches if entry['id'] == transcript][0]            

        print(matches)

        start = start_region[transcript]
        end = end_region[transcript]

        print(start_end)

        if start_end == "start":
            location = f"{start['region_type']} {start['region_number']} "
            region_type = start['region_type']
            region_number = start['region_number']
            other_region_number = end['region_number']
            other_region_type = end['region_type']
            other_location = f"{end['region_type']} {end['region_number']} "
            
            introns = [start['region_number'], end['region_number']]
            position = matches['exons'][start['start-index']] if start['region_type'] == "intron" else matches['exons'][start['start-index'] + 1]
            if location == "intergenic -1 ":
                position = matches['exons'][end['start-index']]

        elif start_end == "end":
            location = f"{end['region_type']} {end['region_number']} "
            region_type = end['region_type']
            region_number = end['region_number']
            other_region_number = start['region_number']
            other_region_type = start['region_type']
            other_location = f"{start['region_type']} {start['region_number']} "
            position = matches['exons'][end['end-index']] if end['region_type'] == "intron" else matches['exons'][end['end-index'] - 1]
            introns = [end['region_number'], start['region_number']]
            if location == "intergenic -1 ":
                position = matches['exons'][start['start-index'] + 1]
            
        
        print(location)
        print(position)
        print(other_region_number)
        print(introns)

        if location == "intergenic -1 ":
            print("INTERGENIC CRYPTIC")
            supp_event = ""
            supp_event_type = ""
            intron = "NA"
            location = "intergenic "
        elif location not in [f"intron {other_region_number} ", f"exon {other_region_number} ", f"exon {other_region_number+1} "]:
            print("but event spans multiple introns")
            if region_type == "exon" and region_number > other_region_number:
                print("region_type: exon")
                introns[0] = introns[0] - 1
            print(introns)
            supp = self._annotate_supplementary(introns)
            supp_event = supp['event']
            supp_event_type = supp["supp_event_type"]
            intron = supp['introns']
        else:
            supp_event = ""
            supp_event_type = ""
            intron = min(introns)

        if start_end == "start":
            if strand == "+":
                event = "donor"
                modifier = 1 if location == f"exon {start['region_number']} " else 0
                distance = query_start - position - modifier
                print(f"The calculation is {query_start} - {position} + {modifier} = {distance}")
            elif strand == "-":
                event = "acceptor"
                modifier = 1 if location in [f"exon {start['region_number']} ", "intergenic "] else 0
                distance = position - query_start + modifier
                print(f"The calculation is {position} - {query_start} + {modifier} = {distance}")
        elif start_end == "end":
            if strand == "+":
                event = "acceptor"
                modifier = 1 if location == f"intron {end['region_number']} " else 0
                distance = query_end - position - modifier
                print(f"The calculation is {query_end} - {position} - {modifier} = {distance}")
            elif strand == "-":
                event = "donor"
                modifier = 1 if location == f"intron {end['region_number']} " else 0
                distance = position - query_end + modifier
                print(f"The calculation is {position} - {query_end}  + {modifier} = {distance}")
        print(event)

        direction = " @ +" if distance > 0 else " @ "

        # exonic acceptor
        # exonic donor

        return {'event': event,
                'cryptic': cryptic,
                'supplementary_event': supp_event,
                'supp_event_type': supp_event_type,
                'location': location,
                'distance': distance,
                'direction': direction,
                'alternate': "",
                'event_type': event,
                'introns': intron}


    def _annotate_unannotated(self, start_region_all, end_region_all):
        print("no start or end matches for transcript identified")
        
        tx = self.coordinates['transcript']

        print(tx)


        start_match = {k: v for k, v in start_region_all.items() if v['transcript'] == tx}
        end_match = {k: v for k, v in end_region_all.items() if v['transcript'] == tx}

        print(start_match)
        print(end_match)
    
        if not start_match:
            region_1 = "intergenic"
            strand_1 = "NA"
        else:
            region_1 = start_match[tx]['region_type']
            strand_1 = start_match[tx]['strand']
        
        if not end_match:
            region_2 = "intergenic"
            strand_2 = "NA"
        else:
            region_2 = end_match[tx]['region_type']
            strand_2 = end_match[tx]['strand']
        
        print(region_1)
        print(region_2)

        print(start_match)
        print(end_match)
        
        if region_1 == "intron" and region_2 == "intron":
            print("both start and end are within introns")
            strand_check = strand_1 == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""
            return {'event': f"junction{strand_in}",
                    'cryptic': "unannotated ",
                    'supplementary_event': "",
                    'supp_event_type': "",
                    'location': "intronic ",
                    'distance': "NA",
                    'direction': "",
                    'alternate': "",
                    'event_type': f"intron{strand_in}",
                    'introns': "NA"}
        
        elif (region_1 == 'exon' and region_2 == 'intron') or (region_2 == 'exon' and region_1 == 'intron'):
            strand_check = strand_1 == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""
            return {'event': f"junction{strand_in}",
                'cryptic': "unannotated ",
                'supplementary_event': "",
                'supp_event_type': "",
                'location': "intronic/exonic ",
                'distance': "NA",
                'direction': "",
                'alternate': "",
                'event_type': f"intron{strand_in}",
                'introns': "NA"}
            
        elif (region_1 == 'exon' and region_2 == 'intergenic') or (region_2 == 'exon' and region_1 == 'intergenic'):
            strand_check = strand_1 == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""
            return {'event': f"junction{strand_in}",
                'cryptic': "unannotated ",
                'supplementary_event': "",
                'supp_event_type': "",
                'location': "intergenic/exonic ",
                'distance': "NA",
                'direction': "",
                'alternate': "",
                'event_type': f"intron{strand_in}",
                'introns': "NA"}
        
        elif (region_1 == 'intron' and region_2 == 'intergenic') or (region_2 == 'intron' and region_1 == 'intergenic'):
            strand_check = strand_1 == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""
            return {'event': f"junction{strand_in}",
                'cryptic': "unannotated ",
                'supplementary_event': "",
                'supp_event_type': "",
                'location': "intergenic/intronic ",
                'distance': "NA",
                'direction': "",
                'alternate': "",
                'event_type': f"intron{strand_in}",
                'introns': "NA"}
        
        elif region_1 == "exon" and region_2 == "exon":
            strand_check = strand_1 == self.coordinates['strand']
            strand_in = " (opposite strand)" if not strand_check else ""    
            return {'event': f"junction{strand_in}",
                'cryptic': "unannotated ",
                'supplementary_event': "",
                'supp_event_type': "",
                'location': "exonic ",
                'distance': "NA",
                'direction': "",
                'alternate': "",
                'event_type': f"intron{strand_in}",
                'introns': "NA"}
        
        else:
            return {'event': "junction",
                'cryptic': "unannotated ",
                'supplementary_event': "",
                'supp_event_type': "",
                'location': "intergenic ",
                'distance': "NA",
                'direction': "",
                'alternate': "",
                'event_type': "intron",
                'introns': "NA"}


    def _annotate_intron_retention(self, start_intron):
        if start_intron['region_type'] == "intron":
            intron = start_intron['region_number']
            cryptic = "intron "
            event = f"{intron} retention"
            event_type = "intron retention"

            return {'event': event,
                    'cryptic': cryptic,
                    'supplementary_event': "",
                    'supp_event_type': "",
                    'location': "",
                    'distance': "",
                    'direction': "",
                    'alternate': "",
                    'event_type': event_type,
                    'introns': intron}

    def _annotate_canonical(self, start_intron):
        if start_intron['region_type'] == "intron":
            intron = start_intron['region_number']
            cryptic = "canonical "
            event = f"exon {intron}-{intron+1} splicing"
            event_type = "canonical"

            return {'event': event,
                    'cryptic': cryptic,
                    'supplementary_event': "",
                    'supp_event_type': "",
                    'location': "",
                    'distance': "",
                    'direction': "",
                    'alternate': "",
                    'event_type': event_type,
                    'introns': intron}

    def _annotate_exon_skipping(self, start_intron, end_intron):
        if start_intron['region_type'] == "intron":
            introns = [start_intron['region_number'], end_intron['region_number']]
            event = f"exon {'-'.join(str(i) for i in range(min(introns)+1, max(introns)+1))} skipping"
            event_type = "exon skipping"
            introns = f"{', '.join(str(i) for i in range(min(introns), max(introns)+1))}"



            return {'event': event,
                    'cryptic': "",
                    'supplementary_event': "",
                    'supp_event_type': "",
                    'location': "NA",
                    'distance': "NA",
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
    
    def read_refgene(self, dataset):
        if dataset == "refseq":
            input_file = "reference/genes.refGene"
        elif dataset == "ensembl":
            raise ValueError("Ensembl annotations are currently not supported. Please try again with 'refseq'.")
        else:
            raise ValueError("Invalid annotation choice. Please select either 'refseq' or 'ensembl' (currently not supported).")
        
        self.refgene = self.read_genepred(input_file, skip_first_column=True)

    def read_genepred(self, input_file, skip_first_column=False):
        """
        GenePred extension format:
        http://genome.ucsc.edu/FAQ/FAQformat.html#GenePredExt

        Column definitions:
        0. string name;                 "Name of gene (usually transcript_id from GTF)"
        1. string chrom;                "Chromosome name"
        2. char[1] strand;              "+ or - for strand"
        3. uint txStart;                "Transcription start position"
        4. uint txEnd;                  "Transcription end position"
        5. uint cdsStart;               "Coding region start"
        6. uint cdsEnd;                 "Coding region end"
        7. uint exonCount;              "Number of exons"
        8. uint[exonCount] exonStarts;  "Exon start positions"
        9. uint[exonCount] exonEnds;    "Exon end positions"
        10. uint id;                    "Unique identifier"
        11. string name2;               "Alternate name (e.g. gene_id from GTF)"
        """

        dataset = []

        with open(input_file, 'r') as infile:
            for line in infile:
                # Skip comments.
                if line.startswith('#'):
                    continue
                row = line.rstrip('\n').split('\t')
                if skip_first_column:
                    row = row[1:]
                # Skip trailing ,
                exon_starts = list(map(int, row[8].split(',')[:-1]))
                exon_ends = list(map(int, row[9].split(',')[:-1]))
                exons = [coord for pair in zip(exon_starts, exon_ends) for coord in pair]
                        
                data = {
                    'chrom': f"chr{row[1]}",
                    'start': int(row[3]),
                    'end': int(row[4]),
                    'id': row[0],
                    'strand': row[2],
                    'cds_start': int(row[5]),
                    'cds_end': int(row[6]),
                    'gene_name': row[11],
                    'exons': exons,
                    'mane': row[15]
                }
                dataset.append(data)

        return dataset