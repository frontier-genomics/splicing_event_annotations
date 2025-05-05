#modified from Himanshu's repo here: https://github.com/kidsneuro-lab/sj_extract/blob/main/cli.py

import argparse
import logging
import os

import src.main as annotator

def main():
    parser = argparse.ArgumentParser(description="Extract split reads for a s/b/cramfile")
    parser.add_argument("-i", "--input", help="File containing split-read coordinates", required=True)
    parser.add_argument("-o", "--output", help="Folder path to write the output file to", required=True)
    parser.add_argument("-g", "--genome", help="Genome build to interpret annotations to: either hg38 or hg19", required=True, default=None)
    parser.add_argument("-d", "--dataset", help="Database to use for annotation (default = refseq)", required=False, default="refseq")
    parser.add_argument("-v", "--verbose", help="Enable verbose logging.", action="store_true")
    # parser.add_argument("-f", "--force", help="Overwrite existing files in output directory.", action="store_true")
    # parser.add_argument("--tempdir", help="Define an alternate directory for temporary intermediate files (default='tempdir/')", required=False, default='tempdir/')

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")
    else:
        logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    logging.info(f"Input:   {args.input}")
    logging.info(f"Genome: {args.genome}")
    logging.info(f"Dataset: {args.dataset}")
    logging.info(f"Output: {args.output}")

    # if not args.output.endswith('/'):
    #     args.output += '/'
    #     logging.debug(f"Added trailing '/' to output folder: {args.output}")
    
    # if not args.tempdir.endswith('/'):
    #     args.tempdir += '/'
    #     logging.debug(f"Added trailing '/' to tempdir folder: {args.output}")

    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file '{args.input}' not found.")
    
    # if not os.path.exists(args.output):
    #     os.makedirs(args.output)
    #     logging.info(f"Created output directory: {args.output}")

    logging.info(f"Annotating splice-junctions from {args.input}")
    annotator.run_annotate(args.input, args.output, args.dataset, args.genome)
    logging.info(f"Completed annotation. Final data saved in {args.output}")
    

if __name__ == "__main__":
    main()