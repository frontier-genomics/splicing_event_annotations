from behave import given, when, then
from src.eventAnnotateProcessorList import EventAnnotateList

@given(u'the splicing events')
def step_impl(context):
    context.splicing_events = []

    for row in context.table.rows:
        row_dict = {
            'chrom': row[1],
            'start': row[2],
            'end': row[3],
            'strand': row[4],
            'type': row[5],
            'transcript': row[6]
        }
        
        context.splicing_events.append(row_dict)


@given(u'the annotation dataset is {dataset}')
def step_impl(context, dataset):
    context.dataset = dataset


@when(u'the list of splicing events is annotated')
def step_impl(context):
    context.annotation_list = EventAnnotateList(
        inputs=context.splicing_events,
        dataset=context.dataset
    )

    context.annotation_list.process()


@then(u'the resulting annotations should be')
def step_impl(context):
    context.annotations = []

    for row in context.table.rows:
        row_dict = {
            'event': row[4],
            'event_type': row[1],
            'introns': row[5],
            'location': row[2],
            'distance_from_authentic': row[3]
        }
        context.annotations.append(row_dict)
    
    print(f"These are the output annotations: \n\n {context.annotation_list.outputs}")
    print(f"\n\n and this is what was expected: \n\n {context.annotations}")

    assert context.annotation_list.outputs == context.annotations