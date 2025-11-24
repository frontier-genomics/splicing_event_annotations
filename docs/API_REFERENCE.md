# API Reference

This document provides detailed API documentation for the SEAR (Splice Event Annotation Repository) library.

## Table of Contents

1. [Main Module](#main-module)
2. [EventAnnotate](#eventannotate)
3. [EventAnnotateList](#eventannotatelist)
4. [TsvAnnotate](#tsvannotate)
5. [BedExpand](#bedexpand)
6. [Data Models](#data-models)

## Main Module

### `splicing_event_annotator.main`

Entry point module for the SEAR library.

#### Functions

##### `run_workflow(input_data, dataset='refseq', genome='hg38', tsv=False, columns=None, output_file=None)`

Main workflow function to annotate splicing events.

**Parameters:**
- `input_data` (dict | list | str): 
  - dict: Single event with keys `chrom`, `start`, `end`, `strand`, `transcript`, `type`
  - list: Multiple events as list of dicts
  - str: Path to TSV file (when `tsv=True`)
- `dataset` (str, optional): Annotation dataset. Options: `'refseq'`, `'ensembl'`. Default: `'refseq'`
- `genome` (str, optional): Genome assembly. Default: `'hg38'`
- `tsv` (bool, optional): Whether input is a TSV file path. Default: `False`
- `columns` (list[int], optional): Column indices for TSV parsing. Default: `[0,1,2,3,4,5]`
- `output_file` (str | Path, optional): Path to write output TSV. Default: `None`

**Returns:**
- `dict | list[dict]`: Annotation results

**Example:**
```python
from splicing_event_annotator.main import run_workflow

# Single event
event = {
    'chrom': 'chr9',
    'start': 34646787,
    'end': 34647088,
    'strand': '+',
    'transcript': 'NM_000155.4',
    'type': 'ir'
}
result = run_workflow(event)
print(result)
# Output: {
#     'event': 'intron 1 retention',
#     'event_type': 'intron retention',
#     'introns': '1',
#     'location': 'NA',
#     'distance_from_authentic': 'NA',
#     'transcript': 'NM_000155.4'
# }

# Multiple events
events = [event1, event2, event3]
results = run_workflow(events)

# From TSV file
results = run_workflow(
    'input.tsv',
    tsv=True,
    output_file='output.tsv'
)
```

---

##### `write_output(output, output_file)`

Write annotation results to a TSV file.

**Parameters:**
- `output` (Iterable[dict]): Annotation results
- `output_file` (str | Path): Output file path

**Raises:**
- `ValueError`: If no annotations were generated

**Example:**
```python
from splicing_event_annotator.main import write_output

results = [annotation1, annotation2, annotation3]
write_output(results, 'annotations.tsv')
```

---

##### `parse_args(argv=None)`

Parse command-line arguments.

**Parameters:**
- `argv` (Sequence[str], optional): Command-line arguments. Default: `sys.argv[1:]`

**Returns:**
- `argparse.Namespace`: Parsed arguments with additional `events` and `output_path` attributes

**Example:**
```python
from splicing_event_annotator.main import parse_args

args = parse_args(['input.tsv', 'output.tsv', '--dataset', 'refseq'])
print(args.dataset)  # 'refseq'
print(args.genome)   # 'hg38' (default)
```

---

##### `main(argv=None)`

CLI entry point.

**Parameters:**
- `argv` (Sequence[str], optional): Command-line arguments

**Example:**
```bash
python -m splicing_event_annotator.main input.tsv output.tsv --dataset refseq --genome hg38
```

---

## EventAnnotate

### `splicing_event_annotator.eventAnnotateProcessor.EventAnnotate`

Core class for annotating individual splicing events.

#### Constructor

##### `EventAnnotate(chrom, start, end, strand, transcript, input_type)`

Initialize an event annotator.

**Parameters:**
- `chrom` (str): Chromosome (e.g., 'chr1', 'chrX')
- `start` (int): Start coordinate (0-based)
- `end` (int): End coordinate (0-based, exclusive)
- `strand` (str): Strand ('+', '-', or '*' for unknown)
- `transcript` (str): Transcript ID (e.g., 'NM_001234.5') or 'NA' for auto-detection
- `input_type` (str): Event type ('sj' for splice junction, 'ir' for intron retention)

**Example:**
```python
from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate

processor = EventAnnotate(
    chrom='chr9',
    start=34646787,
    end=34647088,
    strand='+',
    transcript='NM_000155.4',
    input_type='ir'
)
```

---

#### Methods

##### `process(dataset, genome, get_annotations=True)`

Process the event and generate annotations.

**Parameters:**
- `dataset` (str): Annotation dataset ('refseq' or 'ensembl')
- `genome` (str): Genome assembly (e.g., 'hg38', 'hg19')
- `get_annotations` (bool | list, optional): 
  - `True`: Load annotations from file
  - `False`: Skip annotation loading
  - list: Use pre-loaded annotations. Default: `True`

**Returns:**
- `dict`: Annotation results with keys:
  - `event` (str): Full event description
  - `event_type` (str): Event classification
  - `introns` (str): Affected intron numbers
  - `location` (str): Genomic location
  - `distance_from_authentic` (str): Distance to canonical site (for cryptic events)
  - `transcript` (str): Transcript ID

**Example:**
```python
processor = EventAnnotate('chr9', 34646787, 34647088, '+', 'NM_000155.4', 'ir')
result = processor.process('refseq', 'hg38')
print(result['event'])  # 'intron 1 retention'
```

---

##### `get_mane_transcript(base=1)`

Identify MANE transcript overlapping the event coordinates.

**Parameters:**
- `base` (int, optional): Coordinate system base. Default: `1`

**Returns:**
- `dict`: With keys:
  - `transcript` (str): MANE transcript ID or 'unknown'
  - `gene` (str): Gene symbol or 'unknown'
  - `warning` (str): Warning message or 'none'

**Example:**
```python
processor = EventAnnotate('chr9', 34646787, 34647088, '+', 'NA', 'ir')
mane_info = processor.get_mane_transcript()
print(mane_info['transcript'])  # 'NM_000155.4'
print(mane_info['gene'])        # 'GALT'
```

---

##### `reference_match(start_end)`

Find reference matches for a coordinate.

**Parameters:**
- `start_end` (str): Which coordinate to match ('start' or 'end')

**Returns:**
- `dict`: With keys:
  - `matching_rows` (list): Overlapping transcript records
  - `region_information` (dict): Region details for each transcript

**Raises:**
- `ValueError`: If `start_end` is not 'start' or 'end'

**Example:**
```python
processor = EventAnnotate('chr9', 34646787, 34647088, '+', 'NM_000155.4', 'sj')
processor.refgene = EventAnnotate.read_refgene('refseq', 'hg38')
start_info = processor.reference_match('start')
print(start_info['region_information']['NM_000155.4']['region_type'])  # 'intron'
```

---

##### `get_exon_intron(matching_rows, query, base=0)`

Determine the exon/intron region for a coordinate.

**Parameters:**
- `matching_rows` (list): Transcript records from reference
- `query` (int): Query coordinate
- `base` (int, optional): Coordinate base. Default: `0`

**Returns:**
- `dict`: Region information for each transcript with keys:
  - `transcript` (str): Transcript ID
  - `region_type` (str): 'exon' or 'intron'
  - `region_number` (int): Region number (strand-aware)
  - `match` (bool): Whether coordinate exactly matches boundary
  - `start-index` (int): Start index in exon list
  - `end-index` (int): End index in exon list
  - `strand` (str): Transcript strand

**Example:**
```python
# See reference_match example above
```

---

##### `fetch_transcript_annotations(start_matches_all_tx, end_matches_all_tx)`

Main classification logic for the event.

**Parameters:**
- `start_matches_all_tx` (dict): Start coordinate matching information
- `end_matches_all_tx` (dict): End coordinate matching information

**Returns:**
- `dict`: Complete annotation results

**Example:**
```python
# Called internally by process()
```

---

#### Static Methods

##### `read_refgene(dataset, genome)`

Load reference annotations.

**Parameters:**
- `dataset` (str): 'refseq' or 'ensembl'
- `genome` (str): Genome assembly (e.g., 'hg38')

**Returns:**
- `list[dict]`: Reference annotation records

**Raises:**
- `ValueError`: If dataset is not supported

**Example:**
```python
from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate

annotations = EventAnnotate.read_refgene('refseq', 'hg38')
print(len(annotations))  # Number of transcripts loaded
```

---

##### `read_genepred(input_file, skip_first_column=False)`

Parse a GenePred format file.

**Parameters:**
- `input_file` (str): Path to GenePred file
- `skip_first_column` (bool, optional): Skip first column. Default: `False`

**Returns:**
- `list[dict]`: Parsed records with keys:
  - `chrom` (str): Chromosome
  - `start` (int): Transcription start
  - `end` (int): Transcription end
  - `id` (str): Transcript ID
  - `strand` (str): Strand
  - `cds_start` (int): CDS start
  - `cds_end` (int): CDS end
  - `gene_name` (str): Gene symbol
  - `exons` (list[int]): Flattened exon coordinates
  - `mane` (str): MANE status ('Y' or 'N')

**Example:**
```python
# Using package resources (recommended)
records = EventAnnotate.read_refgene('refseq', 'hg38')

# Or for custom files, use read_genepred directly
from pathlib import Path
records = EventAnnotate.read_genepred(Path('custom_annotations.refGene'), skip_first_column=True)
```

---

## EventAnnotateList

### `splicing_event_annotator.eventAnnotateProcessorList.EventAnnotateList`

Batch processor for multiple splicing events.

#### Constructor

##### `EventAnnotateList(inputs, dataset, genome)`

Initialize batch processor.

**Parameters:**
- `inputs` (list[dict]): List of events (see EventAnnotate for dict format)
- `dataset` (str): Annotation dataset
- `genome` (str): Genome assembly

**Example:**
```python
from splicing_event_annotator.eventAnnotateProcessorList import EventAnnotateList

events = [
    {'chrom': 'chr9', 'start': 34646787, 'end': 34647088, 'strand': '+', 'transcript': 'NM_000155.4', 'type': 'ir'},
    {'chrom': 'chr9', 'start': 34647259, 'end': 34647831, 'strand': '+', 'transcript': 'NM_000155.4', 'type': 'sj'},
]
processor = EventAnnotateList(events, 'refseq', 'hg38')
```

---

#### Methods

##### `load_annotations(dataset, genome)`

Load reference annotations (called automatically by constructor).

**Parameters:**
- `dataset` (str): Annotation dataset
- `genome` (str): Genome assembly

**Returns:**
- `list[dict]`: Reference annotations

---

##### `process()`

Process all events.

**Returns:**
- `list[dict]`: List of annotation results

**Example:**
```python
processor = EventAnnotateList(events, 'refseq', 'hg38')
results = processor.process()
for result in results:
    print(f"{result['event']} - {result['event_type']}")
```

---

## TsvAnnotate

### `splicing_event_annotator.tsvAnnotateProcessor.TsvAnnotate`

TSV file parser.

#### Constructor

##### `TsvAnnotate(tsv, dataset, columns=[0,1,2,3,4,5], header=True)`

Initialize TSV parser.

**Parameters:**
- `tsv` (str): Path to TSV file
- `dataset` (str): Annotation dataset (stored but not used in parsing)
- `columns` (list[int], optional): Column indices for [chrom, start, end, strand, type, transcript]. Default: `[0,1,2,3,4,5]`
- `header` (bool, optional): Whether file has header row. Default: `True`

**Attributes:**
- `tsv` (list[dict]): Parsed events

**Example:**
```python
from splicing_event_annotator.tsvAnnotateProcessor import TsvAnnotate

parser = TsvAnnotate('input.tsv', 'refseq', columns=[0,1,2,3,4,5], header=True)
events = parser.tsv
print(len(events))  # Number of events parsed
```

---

#### Methods

##### `load_tsv(tsv, columns=[0,1,2,3,4,5], header=True)`

Parse TSV file (called automatically by constructor).

**Parameters:**
- `tsv` (str): File path
- `columns` (list[int], optional): Column indices
- `header` (bool, optional): Header presence

**Returns:**
- `list[dict]`: Parsed events

---

## BedExpand

### `splicing_event_annotator.bedExpandProcessor.BedExpand`

BED file parser for curated intron files.

#### Constructor

##### `BedExpand(rows, columns)`

Initialize BED parser.

**Parameters:**
- `rows` (list | pd.DataFrame): BED records
- `columns` (list[str]): Column names

**Example:**
```python
from splicing_event_annotator.bedExpandProcessor import BedExpand

bed_data = [
    ['chr1', 1000, 2000, 'NM_001234_gene_intron_0', 100, '+'],
]
parser = BedExpand(bed_data, ['chr', 'start', 'end', 'name', 'score', 'strand'])
```

---

#### Methods

##### `parse_bed_file()`

Parse BED records.

**Returns:**
- `pd.DataFrame`: Parsed data with columns:
  - `chr` (str): Chromosome
  - `start` (int): Start coordinate
  - `end` (int): End coordinate
  - `transcript` (str): Transcript ID
  - `intron` (int): Intron number (1-based)
  - `strand` (str): Strand

**Example:**
```python
parser = BedExpand(bed_data, columns)
df = parser.parse_bed_file()
print(df[['transcript', 'intron']])
```

---

##### `write_parsed_bed()`

Parse and write to file.

**Returns:**
- `pd.DataFrame`: Parsed data

**Side Effects:**
- Writes to `resources/annotations/curated_introns_sorted.tsv`

---

## Data Models

### SplicingCoordinates

Pydantic model for event coordinates.

**Fields:**
- `chrom` (str): Chromosome in chr[1-X] format
- `start` (int): Start coordinate
- `end` (int): End coordinate
- `strand` (str): Strand ('+', '-', or '*')
- `transcript` (str): Transcript ID
- `type` (str): Event type ('ir' or 'sj')

**Example:**
```python
from splicing_event_annotator.eventAnnotateProcessor import SplicingCoordinates

coords = SplicingCoordinates(
    chrom='chr1',
    start=1000,
    end=2000,
    strand='+',
    transcript='NM_001234.5',
    type='sj'
)
```

---

### AnnotationOutput

Pydantic model for annotation results.

**Fields:**
- `event` (str): Full event description
- `event_type` (str): Event classification
- `introns` (str): Affected intron numbers
- `location` (str): Genomic location or 'NA'
- `distance_from_authentic` (str): Distance value or 'NA'
- `transcript` (str): Transcript ID

**Example:**
```python
from splicing_event_annotator.eventAnnotateProcessor import AnnotationOutput

output = AnnotationOutput(
    event='canonical exon 1-2 splicing',
    event_type='canonical',
    introns='1',
    location='NA',
    distance_from_authentic='NA',
    transcript='NM_001234.5'
)
```

---

## Event Type Classifications

### Output Event Types

| Event Type | Description | Example |
|------------|-------------|---------|
| canonical | Normal exon-exon junction | `canonical exon 1-2 splicing` |
| exon skipping | One or more exons excluded | `exon 2-3 skipping` |
| intron retention | Intron retained in transcript | `intron 1 retention` |
| cryptic acceptor | Novel 3' splice site | `cryptic intron 1 acceptor @ +25` |
| cryptic donor | Novel 5' splice site | `cryptic exon 2 donor @ -15` |
| exon skipping, cryptic acceptor | Combined event | `exon 2 skipping/cryptic intron 3 acceptor @ +10` |
| exon skipping, cryptic donor | Combined event | `exon 2 skipping/cryptic intron 1 donor @ -5` |
| alternate canonical | Matches alternative transcript | `alternate canonical exon 1-2 splicing (NM_001.3)` |
| alternate exon skipping | Alternate transcript pattern | `alternate exon 2 skipping (NM_002.4)` |
| unannotated intron | Novel junction, intronic | `unannotated intronic junction` |
| unannotated intron (opposite strand) | Wrong strand | `unannotated intronic junction (opposite strand)` |
| unknown | Cannot classify | `unknown event (unknown strand)` |

---

## Usage Patterns

### Pattern 1: Single Event Annotation

```python
from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate

# Create processor
processor = EventAnnotate(
    chrom='chr9',
    start=34646787,
    end=34647088,
    strand='+',
    transcript='NM_000155.4',
    input_type='ir'
)

# Get annotation
result = processor.process('refseq', 'hg38')
print(f"Event: {result['event']}")
print(f"Type: {result['event_type']}")
```

---

### Pattern 2: Batch Processing

```python
from splicing_event_annotator.eventAnnotateProcessorList import EventAnnotateList

# Prepare events
events = [
    {'chrom': 'chr9', 'start': 34646787, 'end': 34647088, 'strand': '+', 'transcript': 'NM_000155.4', 'type': 'ir'},
    {'chrom': 'chr9', 'start': 34647259, 'end': 34647831, 'strand': '+', 'transcript': 'NM_000155.4', 'type': 'sj'},
]

# Process batch
processor = EventAnnotateList(events, 'refseq', 'hg38')
results = processor.process()

# Write results
from splicing_event_annotator.main import write_output
write_output(results, 'output.tsv')
```

---

### Pattern 3: TSV Workflow

```python
from splicing_event_annotator.main import run_workflow

# Process TSV file
results = run_workflow(
    input_data='input.tsv',
    dataset='refseq',
    genome='hg38',
    tsv=True,
    output_file='output.tsv'
)
```

---

### Pattern 4: Auto-detect Transcript

```python
from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate

# Use 'NA' to auto-detect MANE transcript
processor = EventAnnotate(
    chrom='chr9',
    start=34646787,
    end=34647088,
    strand='+',
    transcript='NA',  # Auto-detect
    input_type='ir'
)

result = processor.process('refseq', 'hg38')
print(f"Detected transcript: {result['transcript']}")
```

---

## Error Handling

### Common Errors

#### ValueError: Invalid annotation choice

```python
# Wrong dataset name
result = processor.process('invalid', 'hg38')
# Raises: ValueError("Invalid annotation choice. Please select either 'refseq' or 'ensembl' (currently not supported).")
```

**Solution**: Use 'refseq' for dataset parameter.

---

#### ValueError: Missing required columns

```python
# TSV missing columns
args = parse_args(['incomplete.tsv', 'output.tsv'])
# Raises: ValueError("Input file missing required columns: ...")
```

**Solution**: Ensure TSV has columns: chromosome, start, end, strand, event_type, transcript

---

#### ValueError: No annotations were generated

```python
write_output([], 'output.tsv')
# Raises: ValueError("No annotations were generated for the provided input")
```

**Solution**: Check input data and ensure valid events.

---

## Constants

### REQUIRED_COLUMNS

Required columns in input TSV files:

```python
REQUIRED_COLUMNS = (
    "chromosome",
    "start",
    "end",
    "strand",
    "transcript",
    "event_type",
)
```

---

## Best Practices

1. **Use batch processing for multiple events** to avoid reloading annotations
2. **Specify transcripts when known** to avoid MANE lookup overhead
3. **Check event_type in results** to understand classification
4. **Handle 'NA' values** in location and distance fields
5. **Log at INFO level** for production, DEBUG for development
6. **Pre-download reference data** before processing
7. **Validate coordinates** are 0-based, half-open intervals

---

## Version Information

- **Python**: 3.12+
- **Pydantic**: 2.7.0
- **Pandas**: 2.2.1
- **Numpy**: 1.26.4

For updates and issues, see the [GitHub repository](https://github.com/frontier-genomics/splicing_event_annotations).
