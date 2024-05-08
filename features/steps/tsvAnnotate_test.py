from behave import given, when, then
import pandas as pd
import src.tsvAnnotateProcessor as tap


# Annotate test.tsv

@given(u'a tsv {tsv_path}')
def step_impl(context, tsv_path):
    context.tsv_path = tsv_path

@given(u'the annotation dataset is {dataset}')
def step_impl(context, dataset):
    context.dataset = dataset

@when(u'the tsv is annotated')
def step_impl(context):
    context.tsv = tap.TsvAnnotate(context.tsv_path, input_columns = [2,3,4,6])
    context.annotations = context.tsv.annotate(context.dataset)

@then(u'the resulting annotations should be')
def step_impl(context):
    test_data = pd.DataFrame(context.table.rows, columns=context.table.headings)
    output_data = pd.DataFrame(context.annotations, columns=['event'])
    print(f"test_data: \n {test_data['event']}")
    print(f"SEA output: \n {output_data}")
    print(output_data.compare(test_data['event']))
    assert output_data == test_data['event'], f"Expected {test_data['event']}, but got: {output_data}"