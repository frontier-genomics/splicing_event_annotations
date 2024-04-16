# Splice Event Annotation Repository (SEAR)

SEAR is a python library designed to annotate splicing events using data from RefSeq and Ensembl databases. This library provides functionalities to identify and annotate various types of splicing events such as exon skipping, intron retention, alternative splicing, cryptic splice sites, and complex splicing changes.

## What's happening now?
Creating test data by:
1. Take a gene with alternative spliced transcripts
2. Take events in this gene from RNA-sequencing
3. Create annotations for each event by each transcript of the chosen gene
4. Repeat with a few genes to cover different strands
5. Use X chromosome genes to allow subsetting & keep test files small

## Features

- Annotate splicing events using RefSeq and Ensembl databases.
- Identify and classify different types of splicing events.
- Provide a repository of all known human splicing events with annotations
- Support for Python 3.11.6

## Documentation

For detailed documentation, refer to the [official documentation](https://spliceannotator.readthedocs.io/).

## Acknowledgments

- This library utilizes data from RefSeq and Ensembl databases. We acknowledge their contribution to genomic research.
- Special thanks to the developers and contributors of the libraries and tools used in this project.