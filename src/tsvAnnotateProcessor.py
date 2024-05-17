import logging

class TsvAnnotate:

    def __init__(self, tsv, dataset, columns = [0,1,2,3,4,5], header = True):
        self.tsv = self.load_tsv(tsv, columns, header)
        self.dataset = dataset

    def load_tsv(self, tsv, columns = [0,1,2,3,4,5], header = True):
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

        with open(tsv, 'r') as infile:
            if header:
                next(infile)  # Skip the first line
            for line in infile:
                row = line.rstrip('\n').split('\t')
                        
                data = {
                    'chrom': row[columns[0]],
                    'start': int(row[columns[1]]),
                    'end': int(row[columns[2]]),
                    'strand': row[columns[3]],
                    'type': row[columns[4]].lower(),
                    'transcript': row[columns[5]] if len(columns) > 5 else "NA"
                }

                dataset.append(data)

        return dataset