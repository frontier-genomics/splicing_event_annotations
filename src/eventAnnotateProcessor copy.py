import time
import pandas as pd

class EventAnnotate:
    def __init__(self, merged_sjir, merged_sjir_cols, gene, gene_cols):
        self.merged_sjir = pd.DataFrame(merged_sjir, columns=merged_sjir_cols)
        print("self.merged_sjir")
        print(self.merged_sjir)
        self.gene = pd.DataFrame(gene, columns = gene_cols)
        print("self.gene")
        print(self.gene)
        
    def split_name_column(self):
        self.merged_sjir['transcript'] = self.merged_sjir['name'].str.split('_', expand=True).iloc[:, [0, 1]].apply(lambda x: '_'.join(x), axis=1)
        self.merged_sjir['intron'] = self.merged_sjir['name'].str.split('_', expand=True).iloc[:, 3].astype(int) + 1
        self.merged_sjir.drop('name', axis=1, inplace=True)
        print("Intron Retention and Splice Junctions Merged")
        print(self.merged_sjir)
        
    def get_gene_names(self):
        print(self.gene)
        self.merged_sjir_gene = pd.merge(self.merged_sjir, self.gene, left_on=['transcript'], right_on=['tx'], how='inner')
        self.merged_sjir_gene.drop('tx', axis=1, inplace=True)
        print("Genes Annotated")
        print(self.merged_sjir_gene)
        
    def calculate_prop(self):
        self.merged_sjir_gene.sort_values(by=['transcript', 'intron'], inplace=True)
        self.merged_sjir_gene['prop'] = self.merged_sjir_gene.groupby(['transcript', 'intron'])['count'].transform(lambda x: x / x.sum())
        print("Proportions Calculated")
        print(self.merged_sjir_gene)
        
    def event_annotation(self):
        print("Event Annotation")
        self.merged_sjir_gene['result'] = ''
        for index, row in self.merged_sjir_gene.iterrows():
            print(row[['sj_start','ref_start','sj_end','ref_end','event']])
            sj_start = row['sj_start']
            sj_end = row['sj_end']
            ref_start = row['ref_start']
            ref_end = row['ref_end']
            transcript = row['transcript']
            gene = row['gene_name']
            event = row['event']
            sj_intron = row['intron']
            start_match = sj_start == ref_start
            end_match = sj_end == ref_end
            
            if start_match & end_match:
                print("both match")
                print(f"Annotated exon {sj_intron}-{sj_intron+1} splicing or intron {sj_intron} retention")
                if event == 'intron_retention':
                    result = f"Intron {sj_intron} Retention"
                elif event == 'canonical':
                    result = f"Canonical exon {sj_intron}-{sj_intron+1} Splicing"
                print(result)
                self.merged_sjir_gene.loc[index, 'result'] = result
                continue
            
            # Check if non-matching start/end matches another start/end in the same transcript
            if end_match == False:
                print("starts match")
                same_transcript_matches = self.merged_sjir_gene[(self.merged_sjir_gene['transcript'] == transcript) & (self.merged_sjir_gene['ref_end'] == sj_end)]
                #print(same_transcript_matches)
                if not same_transcript_matches.empty:
                    intron = same_transcript_matches.iloc[0]['intron']
                    print(f"Matched another start/end in the same transcript. Intron: {intron}")
                    result = "Exon skipping"
                    print(result)
                    self.merged_sjir_gene.loc[index, 'result'] = result
                    continue
                else:
                    different_transcript_matches = self.merged_sjir_gene[(self.merged_sjir_gene['gene_name'] == gene) & (self.merged_sjir_gene['ref_end'] == sj_end)]
                    if not different_transcript_matches.empty:
                        transcript = different_transcript_matches.iloc[0]['transcript']
                        intron = different_transcript_matches.iloc[0]['intron']
                        print(f"Matched another start/end in a different transcript. Transcript: {transcript}, Intron: {intron}")
                        result = "Exon skipping to alternate exon"
                        print(result)
                        self.merged_sjir_gene.loc[index, 'result'] = result
                        continue
                    else:
                        if sj_end < ref_end:
                            print(f"Positioned before the annotated intron end: Intronic Cryptic Acceptor -{ref_end - sj_end}")
                            result = f"Intronic Cryptic Acceptor -{ref_end - sj_end}"
                            print(result)
                            self.merged_sjir_gene.loc[index, 'result'] = result
                            continue
                        else:
                            print(f"Positioned after the annotated intron end: Exonic/Intronic Cryptic Acceptor +{sj_end - ref_end}")
                            result = f"Exonic/Intronic Cryptic Acceptor +{ref_end - sj_end}"
                            print(result)
                            self.merged_sjir_gene.loc[index, 'result'] = result
                            continue
                            #introns = self.merged_sjir_gene[(self.merged_sjir_gene['transcript'] == row['transcript']) & (self.merged_sjir_gene['intron'] < sj_start)]['intron']
                            #if not introns.empty:
                            #    ref_start_intron = introns.max()
                            #    print(f"Positioned between the annotated start/end. Difference: {ref_start - sj_start}, Intron: {ref_start_intron}")
                            #    continue
                
            elif start_match == False:
                print("ends match")
                same_transcript_matches = self.merged_sjir_gene[(self.merged_sjir_gene['transcript'] == transcript) & (self.merged_sjir_gene['ref_start'] == sj_start)]
                #print(same_transcript_matches)
                if not same_transcript_matches.empty:
                    intron = same_transcript_matches.iloc[0]['intron']
                    print(f"Matched another start/end in the same transcript. Intron: {intron}")
                    continue
                else:
                    different_transcript_matches = self.merged_sjir_gene[(self.merged_sjir_gene['gene_name'] == gene) & (self.merged_sjir_gene['ref_start'] == sj_start)]
                    if not different_transcript_matches.empty:
                        transcript = different_transcript_matches.iloc[0]['transcript']
                        intron = different_transcript_matches.iloc[0]['intron']
                        print(f"Matched another start/end in a different transcript. Transcript: {transcript}, Intron: {intron}")
                        continue
                    else:
                        if sj_start > ref_start:
                            print(f"Positioned after the annotated intron start: Intronic Cryptic Donor +{sj_start - ref_start}")
                            result = f"Intronic Cryptic Donor +{sj_start - ref_start}"
                            print(result)
                            self.merged_sjir_gene.loc[index, 'result'] = result
                            continue
                        else:
                            print(f"Positioned before the annotated intron start: Exonic/Intronic Cryptic Donor -{ref_start - sj_start}")
                            result = f"Exonic/Intronic Cryptic Donor -{ref_start - sj_start}"
                            print(result)
                            self.merged_sjir_gene.loc[index, 'result'] = result
                            continue
        print(self.merged_sjir_gene)
                            #introns = self.merged_sjir_gene[(self.merged_sjir_gene['transcript'] == row['transcript']) & (self.merged_sjir_gene['intron'] < sj_start)]['intron']
                            #if not introns.empty:
                            #    ref_start_intron = introns.max()
                            #    print(f"Positioned between the annotated start/end. Difference: {ref_start - sj_start}, Intron: {ref_start_intron}")
                            #    continue
            
            # Check if non-matching start/end matches another start/end in a different transcript
            # different_transcript_matches = self.merged_sjir_gene[(self.merged_sjir_gene['gene_name'] == gene) & ((self.merged_sjir_gene['sj_start'] == sj_start) | (self.merged_sjir_gene['sj_end'] == sj_end))]
            # if not different_transcript_matches.empty:
            #     transcript = different_transcript_matches.iloc[0]['transcript']
            #     intron = different_transcript_matches.iloc[0]['intron']
            #     print(f"Matched another start/end in a different transcript. Transcript: {transcript}, Intron: {intron}")
            #     continue
            
            # Check if non-matching start/end is positioned between the annotated start/end
            # if sj_start != ref_start:
            #     if sj_start > ref_start:
            #         print(f"Positioned between the annotated start/end. Difference: {ref_start - sj_start}")
            #         continue
            #     else:
            #         introns = self.merged_sjir_gene[(self.merged_sjir_gene['transcript'] == row['transcript']) & (self.merged_sjir_gene['intron'] < sj_start)]['intron']
            #         if not introns.empty:
            #             ref_start_intron = introns.max()
            #             print(f"Positioned between the annotated start/end. Difference: {ref_start - sj_start}, Intron: {ref_start_intron}")
            #             continue
            
            #if sj_end != ref_end:
            #    if sj_end < ref_end:
            #        print(f"Positioned between the annotated start/end. Difference: {ref_end - sj_end}")
            #        continue
            #    else:
            #        introns = self.merged_sjir_gene[(self.merged_sjir_gene['transcript'] == row['transcript']) & (self.merged_sjir_gene['intron'] > sj_end)]['intron']
            #        if not introns.empty:
            #            ref_end_intron = introns.min()
            #            print(f"Positioned between the annotated start/end. Difference: {ref_end - sj_end}, Intron: {ref_end_intron}")
            #            continue
    
    #def refactor_output(self):

    
    def run(self):
        # get the start time
        start_time = time.time()        
        # get the end time
        end_time = time.time()
        # get the execution time
        self.runtime = end_time - start_time
        print('Total execution time:', self.runtime, 'seconds')