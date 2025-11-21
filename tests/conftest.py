import pytest
import sys
from pathlib import Path

# Add the parent directory to sys.path to allow imports
sys.path.insert(0, str(Path(__file__).parent.parent))


@pytest.fixture(scope="session")
def refgene_annotations():
    """Load RefSeq annotations once for all tests.
    
    This session-scoped fixture loads the annotation file once and caches it
    for reuse across all tests, significantly improving test performance.
    """
    from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate
    
    # Load the annotations directly using the static method (no instantiation needed)
    annotations = EventAnnotate.read_refgene('refseq', 'hg38')
    
    return annotations


