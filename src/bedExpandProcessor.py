import pandas as pd

class BedExpand:
    def __init__(self, rows, columns):
        if isinstance(rows, pd.DataFrame):
            self.bed_file = rows.copy()
        else:
            self.bed_file = pd.DataFrame(rows, columns=columns)
        
    def parse_bed_file(self):
        df = self.bed_file
        df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand']
        split_names = df['name'].str.split('_')
        df['transcript'] = split_names.str[0] + "_" + split_names.str[1]
        df['intron'] = split_names.str[3].astype(int) + 1
        
        return df[['chr', 'start', 'end', 'transcript', 'intron', 'strand']]
    
    def write_parsed_bed(self):
        parsed_bed = self.parse_bed_file()
        parsed_bed.to_csv("resources/annotations/curated_introns_sorted.tsv", sep='\t', header=True, index=False)        
        return parsed_bed
