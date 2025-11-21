from splicing_event_annotator.eventAnnotateProcessorList import EventAnnotateList
from splicing_event_annotator.tsvAnnotateProcessor import TsvAnnotate
import logging

logging.basicConfig(level=logging.DEBUG)

tsv_processor = TsvAnnotate("resources/test_input.tsv", "refseq", [2, 3, 4, 6, 25])

input = tsv_processor.tsv

processor = EventAnnotateList(input, "refseq", "hg38")

processor.process()
