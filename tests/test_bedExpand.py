import pytest
from pytest_bdd import given, when, then, scenarios
from splicing_event_annotator.bedExpandProcessor import BedExpand
import pandas as pd

scenarios('bedExpand.feature')

@given('the RefSeq MANE transcript bed file is', target_fixture='bed_file')
def bed_file_input(datatable):
    """Parse the input bed file from the data table."""
    # datatable is provided by pytest-bdd 8.x as list of lists
    # First row is headers, rest are data rows
    headers = datatable[0]
    rows = datatable[1:]
    
    return BedExpand(rows, headers)

@when('the RefSeq MANE transcript bed file is parsed', target_fixture='expanded_result')
def parse_bed_file(bed_file):
    """Parse the bed file."""
    return bed_file.parse_bed_file()

@then('the RefSeq MANE transcript bed file should be')
def verify_bed_file(expanded_result, datatable):
    """Verify the parsed bed file matches expected output."""
    expanded = pd.DataFrame(expanded_result).astype(str)
    print(expanded)
    
    # datatable is list of lists: first row is headers, rest are data
    headers = datatable[0]
    rows = datatable[1:]
    
    test_data = pd.DataFrame(rows, columns=headers).astype(str)
    print(test_data)
    
    difference = expanded.compare(test_data)
    print(difference)
    assert expanded.equals(test_data)

