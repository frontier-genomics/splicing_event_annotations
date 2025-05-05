from behave import given, when, then
import src.main as main
import logging

@given(u'the splicing event tsv {tsv}')
def step_impl(context, tsv):
    context.tsv = tsv

@given(u'the annotation dataset for the tsv is {dataset}')
def step_impl(context, dataset):
    context.dataset = dataset

@when(u'the tsv of splicing events is annotated')
def step_impl(context):
    context.output = main.run_workflow(context.tsv, context.dataset, genome = "hg38", tsv = True, columns = [2,3,4,6,25])
    print(context.output)
    main.write_output(context.output, 'resources/test_output.tsv')

@then(u'the resulting annotated tsv should be {output}')
def step_impl(context, output):
    raise NotImplementedError(u'STEP: Then the resulting annotated tsv should be "resources/test_output.tsv"')