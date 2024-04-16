import pandas as pd
import numpy as np

class EventAnnotate:
    """
    Class representing an event annotation.

    Attributes:
        coordinates (dict): A dictionary containing the coordinates of the event.
        annotation_choice (str): The choice of annotation to use.
        annotation (pd.DataFrame): The annotation data.

    Methods:
        matcher: Matches the event coordinates with the annotation data.
    """

    def __init__(self, chrom, start, end, strand, type, annotation_choice):
        """
        Initializes an EventAnnotate object.

        Args:
            chrom (str): The chromosome of the event.
            start (int): The start position of the event.
            end (int): The end position of the event.
            strand (str): The strand of the event.
            type (str): The type of the event.
            annotation_choice (str): The choice of annotation to use.
        
        Raises:
            ValueError: If the annotation_choice is invalid.
        """
        self.coordinates = {
            'chrom': str(chrom),
            'start': np.int64(start),
            'end': np.int64(end),
            'strand': str(strand),
            'type': str(type)
        }
        self.annotation_choice = annotation_choice

        if self.annotation_choice == "refseq_mane":
            self.annotation = pd.read_csv("resources/annotations/refseq_mane_introns_sorted.bed", sep='\t', header=None)
        elif self.annotation_choice == "refseq_curated":
            self.annotation = pd.read_csv("resources/annotations/curated_introns_sorted.bed", sep='\t', header=None)
        else:
            raise ValueError("Invalid annotation choice. Please select either 'refseq_mane' or 'refseq_curated'.")
        
        self.annotation.columns = ['chrom', 'start', 'end', 'intron', 'score', 'strand']


    def matcher(self, start_end, base = 1):
        """
        Matches the event coordinates with the annotation data.

        Args:
            start_end (str): The choice of start or end coordinate to match.
            base (int, optional): The base to subtract from the start coordinate. Defaults to 1.
        
        Raises:
            ValueError: If the start_end choice is invalid.
        """
        if start_end not in ['start', 'end']:
            raise ValueError("Invalid start_end choice. Please select either 'start' or 'end'.")
        
        if start_end == "start":
            query = self.coordinates[start_end] - base
        elif start_end == "end":
            query = self.coordinates[start_end]
        
        chrom = self.coordinates['chrom']
        
        print(query)
        print(self.annotation[start_end])
        
        matching_rows = self.annotation[
            (self.annotation[start_end] == query) &
            (self.annotation['chrom'] == chrom)
        ]
        
        match = pd.DataFrame(matching_rows)
        print(match)