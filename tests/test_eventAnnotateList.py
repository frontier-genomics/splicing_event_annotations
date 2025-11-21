import pytest
from pytest_bdd import given, when, then, scenarios, parsers
from splicing_event_annotator.eventAnnotateProcessorList import EventAnnotateList
import logging

logging.basicConfig(level=logging.DEBUG)

scenarios('eventAnnotateList.feature')

@pytest.fixture
def splicing_events():
    return []

@pytest.fixture
def dataset_name():
    return {}

@pytest.fixture
def annotation_list_result():
    return None

@given(parsers.parse('the splicing events'))
def given_splicing_events(splicing_events, datatable):
    for row in datatable:
        row_values = list(row.values())
        row_dict = {
            'chrom': row_values[1],
            'start': row_values[2],
            'end': row_values[3],
            'strand': row_values[4],
            'type': row_values[5],
            'transcript': row_values[6]
        }
        splicing_events.append(row_dict)

@given(parsers.parse('the annotation dataset is {dataset}'))
def given_annotation_dataset(dataset_name, dataset):
    dataset_name['value'] = dataset

@when(parsers.parse('the list of splicing events is annotated'), target_fixture='annotation_list_result')
def annotate_splicing_events(splicing_events, dataset_name):
    annotation_list = EventAnnotateList(
        inputs=splicing_events,
        dataset=dataset_name['value'],
        genome='hg38'
    )
    annotation_list.process()
    return annotation_list

@then(parsers.parse('the resulting annotations should be'))
def verify_annotations(annotation_list_result, datatable):
    annotations = []
    for row in datatable:
        row_values = list(row.values())
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
