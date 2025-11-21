import pytest
from pytest_bdd import given, when, then, scenarios, parsers
from splicing_event_annotator.bedExpandProcessor import BedExpand
import pandas as pd

scenarios('bedExpand.feature')

@given(parsers.parse('the RefSeq MANE transcript bed file is'), target_fixture='bed_file')
def bed_file_input(datatable):
    """Parse the input bed file from the data table."""
    rows = []
    for row in datatable:
        rows.append(list(row.values()))
    headings = list(datatable[0].keys())
    return BedExpand(rows, headings)

@when(parsers.parse('the RefSeq MANE transcript bed file is parsed'), target_fixture='expanded_result')
def parse_bed_file(bed_file):
    """Parse the bed file."""
    return bed_file.parse_bed_file()

@then(parsers.parse('the RefSeq MANE transcript bed file should be'))
def verify_bed_file(expanded_result, datatable):
    """Verify the parsed bed file matches expected output."""
    expanded = pd.DataFrame(expanded_result).astype(str)
    print(expanded)
    
    rows = []
    for row in datatable:
        rows.append(list(row.values()))
    headings = list(datatable[0].keys())
    test_data = pd.DataFrame(rows, columns=headings).astype(str)
    print(test_data)
    
    difference = expanded.compare(test_data)
    print(difference)
    assert expanded.equals(test_data)
