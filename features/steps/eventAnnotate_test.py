from behave import given, when, then
from src.eventAnnotateProcessor import EventAnnotate

@given(u'a splicing event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
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
        transcript = context.input['transcript'],
        type = context.input['type'],
        annotation_choice = context.input['annotation']
    )
    context.annotation.get_annotations()
    context.start = context.annotation.matcher('start')
    context.end = context.annotation.matcher('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.create_annotations(context.start, context.end)
    

@then(u'the result should be event {event} for transcript {transcript}')
def step_impl(context, event, transcript):
    print(context.create_annotations['event'])
    print(context.create_annotations['transcript'])
    print(event)
    print(transcript)
    assert context.create_annotations['event'] == event
    assert context.create_annotations['transcript'] == transcript





@given(u'another splicing event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
    
@given(u'the selected transcript annotation is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the splicing events are annotated with the curated annotations')
def step_impl(context):
    context.annotation = EventAnnotate(
        chrom = context.input['chrom'],
        start = context.input['start'],
        end = context.input['end'],
        strand = context.input['strand'],
        transcript = context.input['transcript'],
        type = context.input['type'],
        annotation_choice = context.input['annotation']
    )
    context.annotation.get_annotations()
    context.start = context.annotation.matcher('start')
    context.end = context.annotation.matcher('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.create_annotations(context.start, context.end)
    

@then(u'the resulting annotations should be event {event} for transcript {transcript}')
def step_impl(context, event, transcript):
    print(context.create_annotations['event'])
    print(context.create_annotations['transcript'])
    print(event)
    print(transcript)
    assert context.create_annotations['event'] == event
    assert context.create_annotations['transcript'] == transcript