import concurrent.futures
import pysam
import logging
import os

class samReader:
    def __init__(self, samfile, tempdir, reference_file, threads, chrom_lengths = "resources/chrom_lengths.tsv", force = False):
        self.file = samfile
        self.extension = os.path.splitext(samfile)[1][1:]
        self.first_letter = self.extension[0]
        self.tempdir = f"{tempdir}{os.path.splitext(os.path.basename(samfile))[0]}/"
        self.prefix = os.path.splitext(os.path.basename(samfile))[0]
        self.reference = reference_file
        self.chromosomes = self.get_chrom_lengths(chrom_lengths)
        self.divide_genome_into_intervals()
        self.threads = int(threads)
        self.force = force

        if not os.path.exists(self.tempdir):
            os.makedirs(self.tempdir)
        elif os.path.exists(self.tempdir) and self.force:
            logging.debug(f"Force flag detected. Deleting existing tempdir: {self.tempdir}")
            os.system(f"rm -r {self.tempdir}")
            os.makedirs(self.tempdir)
        else:
            logging.error(f"Output directory '{self.tempdir}' already exists. Use the -f flag to overwrite existing files.")
            raise FileExistsError(f"Output directory '{self.tempdir}' already exists. Use the -f flag to overwrite existing files.")

        logging.debug(f"Detected input to be of file type: '.{self.extension}'")
        if self.extension == 'cram':
            logging.debug(f"Using FASTA reference located at: {self.reference}")
        
        logging.debug(f"Saving intermediate files to: {self.tempdir}")

    def process_interval(self, interval, mode='parallel'):
        self.first_letter = '' if self.first_letter == 's' else self.first_letter
        with pysam.AlignmentFile(self.file, "r" + self.first_letter, reference_filename = self.reference) as samfile:
            if mode != 'parallel':
                output_file_path = f"{self.tempdir}{self.prefix}.tsv"
                logging.debug(f"STARTED PROCESS - Extracting all splice-junctions")
            else:
                output_file_path = f"{self.tempdir}{self.prefix}-{interval[0]}_{interval[1]}_{interval[2]}.tsv"
                logging.debug(f"STARTED PROCESS - Extracting splice-junctions for interval: {interval}")
            with open(output_file_path, "w") as tsv_file:
                tsv_file.write("flag\treference_name\tmapping_quality\tNH\tAS\tjI\tjM\tXS\n")
                for read in samfile.fetch(contig = interval[0], start = interval[1], stop = interval[2]):
                    row_data = [
                        read.flag,
                        read.reference_name,
                        read.mapping_quality,
                        read.get_tag("NH") if read.has_tag("NH") else None,
                        read.get_tag("AS") if read.has_tag("AS") else None,
                        ','.join(str(element) for element in read.get_tag("jI")) if read.has_tag("jI") else None,
                        ','.join(str(element) for element in read.get_tag("jM")) if read.has_tag("jM") else None,
                        read.get_tag("XS") if read.has_tag("XS") else None
                    ]
                    if row_data[5] != '-1':
                        tsv_file.write('\t'.join(str(element) for element in row_data) + '\n')
            if mode != 'parallel':
                logging.debug(f"PROCESS ENDED - Extracting splice-junctions")
            else:
                logging.debug(f"PROCESS ENDED - Extracting splice-junctions for interval: {interval}")
                
        
    def divide_genome_into_intervals(self, interval_size=100000000):
        intervals = []
        for chr_name, chr_length in self.chromosomes.items():
            for start in range(1, chr_length + 1, interval_size):
                end = min(start + interval_size - 1, chr_length)
                intervals.append((chr_name, start, end))
        self.intervals = intervals

    def extract_sj_counts(self):
        if(self.threads == 1):
            mode = 'single'
            self.process_interval(interval=[None, None, None], mode=mode)
        else:
            mode = 'parallel'
            with concurrent.futures.ProcessPoolExecutor(max_workers=self.threads) as executor:
                [executor.submit(self.process_interval, interval, mode) for interval in self.intervals]

    def get_chrom_lengths(self, chrom_lengths_file):
        chrom_lengths = {}
        with open(chrom_lengths_file, 'r') as file:
            for line in file:
                parts = line.strip().split('\t')  # Remove newline and split by tab
                if len(parts) == 2:  # Ensure the line has exactly two parts
                    chrom, length = parts
                    chrom_lengths[chrom] = int(length)  # Convert length to integer and add to dictionary
        return chrom_lengths


'''
import pysam
import os
import src.readSAM as rs
import time
test = rs.samReader("/host/host_mnt/Volumes/research-data/PRJ-GDT/NGS_data/Rhett/alignments/results/160_CHD7_P_PAX/160_CHD7_P_PAX_var1_6bp/160_CHD7_P_PAX_var1_6bp.sorted.cram", "tempdir/test_parallel_cram", "/host/host_mnt/Users/rmar4592/Documents/local_ref/hg38.fa")
start = time.time()
test.extract_sj_counts()
end = time.time()
print(end - start)

python clip.py -i /host/host_mnt/Volumes/research-data/PRJ-GDT/NGS_data/Rhett/alignments/results/160_CHD7_P_PAX/160_CHD7_P_PAX_var1_6bp/160_CHD7_P_PAX_var1_6bp.sorted.cram -r /host/host_mnt/Users/rmar4592/Documents/local_ref/hg38.fa -t 4 -v

datamash --sort --headers --group 6,1 count 2  < temp/head_test.tsv > head_test_gb_jI_flag.tsv
datamash --sort -g 1,2 sum 3 <  head_test_gb_jI_flag.tsv

use R code to do the next part

datamash --sort  --headers --group 6,1 count 2  < temp/test.tsv > temp/test_grouped_ji.tsv
cat single.tsv >> merged.tsv
cat multi.tsv >> merged.tsv
datamash --sort -g 1,2,3 sum 4 < temp/merged.tsv > temp/merged_sum.tsv
awk '{gsub(/,/, "\t", $0); print}' temp/merged_sum.tsv > temp/merged_sum_tab.tsv


test.divide_genome_into_intervals()
test.file
test.extension
test.first_letter
'''