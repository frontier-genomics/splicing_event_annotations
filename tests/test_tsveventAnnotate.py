import pytest
from pytest_bdd import given, when, then, scenarios, parsers
import splicing_event_annotator.main as main
import logging

scenarios('tsveventAnnotate.feature')

@given(parsers.parse('the splicing event tsv {tsv}'), target_fixture='tsv_path')
def given_tsv(tsv):
    """Store the TSV path."""
    return tsv

@given(parsers.parse('the annotation dataset for the tsv is {dataset}'), target_fixture='dataset_info')
def given_dataset(dataset):
    """Store the dataset name."""
    return dataset

@when('the tsv of splicing events is annotated', target_fixture='output_result')
def annotate_tsv(tsv_path, dataset_info):
    """Annotate the TSV file."""
    output = main.run_workflow(tsv_path, dataset_info, tsv=True, columns=[2,3,4,6,25])
    print(output)
    main.write_output(output, 'output.tsv')
    return output

@then(parsers.parse('the resulting annotated tsv should be {output}'))
def verify_tsv_output(output_result, output):
    """Verify the TSV output."""
    raise NotImplementedError(u'STEP: Then the resulting annotated tsv should be "resources/test_output.tsv"')

