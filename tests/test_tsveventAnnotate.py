import pytest
from pytest_bdd import given, when, then, scenarios, parsers
import splicing_event_annotator.main as main
import logging

scenarios('tsveventAnnotate.feature')

@pytest.fixture
def tsv_path():
    return {}

@pytest.fixture
def dataset_info():
    return {}

@pytest.fixture
def output_result():
    return None

@given(parsers.parse('the splicing event tsv {tsv}'))
def given_tsv(tsv_path, tsv):
    tsv_path['value'] = tsv

@given(parsers.parse('the annotation dataset for the tsv is {dataset}'))
def given_dataset(dataset_info, dataset):
    dataset_info['value'] = dataset

@when(parsers.parse('the tsv of splicing events is annotated'), target_fixture='output_result')
def annotate_tsv(tsv_path, dataset_info):
    output = main.run_workflow(tsv_path['value'], dataset_info['value'], tsv=True, columns=[2,3,4,6,25])
    print(output)
    main.write_output(output, 'output.tsv')
    return output

@then(parsers.parse('the resulting annotated tsv should be {output}'))
def verify_tsv_output(output_result, output):
    raise NotImplementedError(u'STEP: Then the resulting annotated tsv should be "resources/test_output.tsv"')
