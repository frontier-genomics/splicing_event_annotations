from behave import given, when, then
from src.eventAnnotateProcessor import EventAnnotate

@given(u'a splicing event with chrom {chrom}, start {start}, end {end}, strand {strand}, and type {type}')
def step_impl(context, chrom, start, end, strand, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'type': type
        }
    
@given(u'the transcript annotation is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the splicing events are annotated')
def step_impl(context):
    context.annotation = EventAnnotate(
        chrom = context.input['chrom'],
        start = context.input['start'],
        end = context.input['end'],
        strand = context.input['strand'],
        type = context.input['type'],
        annotation_choice = context.input['annotation']
    )
    context.start = context.annotation.matcher('start')
    context.end = context.annotation.matcher('end')

@then(u'the result should be event {event} for transcript {transcript}')
def step_impl(context, event, transcript):
    raise NotImplementedError(u'the result should be event {event} for transcript {transcript}')