from src.eventAnnotateProcessor import EventAnnotate
from src.eventAnnotateProcessorList import EventAnnotateList
from src.tsvAnnotateProcessor import TsvAnnotate
import csv


def run_annotate(input_file, output_file, dataset, genome):

    annotation_database = read_refgene(dataset, genome)

    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')
        
        for row in reader:        
            input_data = {'chrom': row[0], 'start': row[1], 'end': row[2], 'strand': row[3], 'transcript': row[4], 'type': row[5]}
            processor = EventAnnotate(input_data['chrom'], input_data['start'], input_data['end'], input_data['strand'], input_data['transcript'], input_data['type'])
            annotations = processor.process(dataset, genome, annotation_database)
            
            writer.writerow(row + list(annotations.values()))

def read_refgene(dataset, genome):
    if dataset == "refseq":
        input_file = f"reference/{genome}/genes.refGene"
    elif dataset == "ensembl":
        raise ValueError("Ensembl annotations are currently not supported. Please try again with 'refseq'.")
    else:
        raise ValueError("Invalid annotation choice. Please select either 'refseq' or 'ensembl' (currently not supported).")
    
    refgene = read_genepred(input_file, skip_first_column=True)
    return refgene

def read_genepred(input_file, skip_first_column=False):
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

def run_workflow(input, dataset, genome, tsv = False, columns = [0,1,2,3,4,5], output_file = ""):

    if tsv == True:
        tsv_processor = TsvAnnotate(input, dataset, columns)
        input = tsv_processor.tsv

    if isinstance(input, list):
        processor = EventAnnotateList(input, dataset, genome)
    else:
        processor = EventAnnotate(input['chrom'], input['start'], input['end'], input['strand'], input['transcript'], input['type'])

    annotations = processor.process()

    if output_file != "":
        write_output(annotations, output_file)

    return annotations


def write_output(output, output_file):
    keys = output[0].keys()

    with open(output_file, 'w', newline='') as file:
        dict_writer = csv.DictWriter(file, keys, delimiter='\t')
        dict_writer.writeheader()
        dict_writer.writerows(output)