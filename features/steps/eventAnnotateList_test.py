from behave import given, when, then
import src.eventAnnotateProcessorList

@given(u'the splicing events')
def step_impl(context):
    pass

@given(u'the annotation dataset is {dataset}')
def step_impl(context, dataset):
    context.dataset = dataset

@when(u'the list of splicing events is annotated')
def step_impl(context):
    pass

@then(u'the resulting annotations should be')
def step_impl(context):
    print(context.table)
    assert context.table == 'test'