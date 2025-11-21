import pytest
from pytest_bdd import given, when, then, scenarios, parsers
from splicing_event_annotator.eventAnnotateProcessorList import EventAnnotateList
import logging

logging.basicConfig(level=logging.DEBUG)

scenarios('eventAnnotateList.feature')

@given('the splicing events')
def given_splicing_events(datatable):
    """Store splicing events from datatable."""
    # datatable is list of lists with headers in first row
    headers = datatable[0]
    events = []
    
    for row_values in datatable[1:]:
        row_dict = {
            'chrom': row_values[1],
            'start': row_values[2],
            'end': row_values[3],
            'strand': row_values[4],
            'type': row_values[5],
            'transcript': row_values[6]
        }
        events.append(row_dict)
    
    return events

@given(parsers.parse('the annotation dataset is {dataset}'), target_fixture='dataset_name')
def given_annotation_dataset(dataset):
    """Store the annotation dataset name."""
    return dataset

@when('the list of splicing events is annotated', target_fixture='annotation_list_result')
def annotate_splicing_events(given_splicing_events, dataset_name):
    """Annotate the list of splicing events."""
    annotation_list = EventAnnotateList(
        inputs=given_splicing_events,
        dataset=dataset_name,
        genome='hg38'
    )
    annotation_list.process()
    return annotation_list

@then('the resulting annotations should be')
def verify_annotations(annotation_list_result, datatable):
    """Verify the annotations match expected output."""
    # datatable is list of lists with headers in first row
    annotations = []
    for row_values in datatable[1:]:
        row_dict = {
            'event': row_values[4],
            'event_type': row_values[1],
            'introns': row_values[5],
            'location': row_values[2],
            'distance_from_authentic': row_values[3]
        }
        annotations.append(row_dict)
    
    print(f"These are the output annotations: \n\n {annotation_list_result.outputs}")
    print(f"\n\n and this is what was expected: \n\n {annotations}")
    
    assert annotation_list_result.outputs == annotations

