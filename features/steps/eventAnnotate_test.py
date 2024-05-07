from behave import given, when, then
from src.eventAnnotateProcessor import EventAnnotate

@given(u'a splicing event with chrom {chrom}, start {start}, end {end}, and strand {strand}')
def step_impl(context, chrom, start, end, strand):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': "",
            'type': ""
        }

@when(u'the overlapping transcripts are fetched')
def step_impl(context):
    context.annotation = EventAnnotate(
        chrom = context.input['chrom'],
        start = context.input['start'],
        end = context.input['end'],
        strand = context.input['strand'],
        transcript = context.input['transcript'],
        type = context.input['type'],
        annotation_choice=""
    )
    context.transcript = context.annotation.get_mane_transcript()


@then(u'the result should be transcript {transcript}')
def step_impl(context, transcript):
    print(context.transcript['transcript'])
    print(transcript)
    assert context.transcript['transcript'] == transcript





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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    

@then(u'the resulting annotations should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"


# CANONICAL SPLICING

@given(u'a canonical splicing event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
    
@given(u'the annotation dataset to be used for annotating canonical splicing is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the canonical splicing events are annotated with the annotation dataset')
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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    

@then(u'the resulting annotations of canonical splicing should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"



# EXON SKIPPING

@given(u'an exon skipping event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
@given(u'the annotation dataset to be used for annotating exon skipping is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the exon skipping events are annotated with the annotation dataset')
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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of exon skipping should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"


# INTRON RETENTION

@given(u'an intron retention event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
@given(u'the annotation dataset to be used for annotating intron retention is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the intron retention events are annotated with the annotation dataset')
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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of intron retention should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"


# CRYPTIC DONOR

@given(u'a cryptic donor event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
@given(u'the annotation dataset to be used for annotating cryptic donors is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the cryptic donor events are annotated with the annotation dataset')
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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of cryptic donors should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"


# CRYPTIC ACCEPTOR

@given(u'a cryptic acceptor event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
@given(u'the annotation dataset to be used for annotating cryptic acceptors is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the cryptic acceptor events are annotated with the annotation dataset')
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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of cryptic acceptors should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"


# SKIPPING WITH A CRYPTIC (SKIP CRYPS)

@given(u'a skip-cryp event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
@given(u'the annotation dataset to be used for annotating skip-cryps is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the skip-cryp events are annotated with the annotation dataset')
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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of skip-cryps should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"


# UNANNOTATED JUNCTIONS

@given(u'an unannotated junction with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
@given(u'the annotation dataset to be used for annotating unannotated junctions is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the unannotated junctions are annotated with the annotation dataset')
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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of unannotated junctions should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"


# ALTERNATE SPLICING

@given(u'an alternate splicing event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}')
def step_impl(context, chrom, start, end, strand, transcript, type):
    context.input = {
            'chrom': chrom,
            'start': start,
            'end': end,
            'strand': strand,
            'transcript': transcript,
            'type': type
        }
@given(u'the annotation dataset to be used for annotating alternate splicing is {annotation}')
def step_impl(context, annotation):
    context.input['annotation'] = annotation

@when(u'the alternate splicing events are annotated with the annotation dataset')
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
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of alternate splicing should be event {event}')
def step_impl(context, event):
    print(str(context.create_annotations))
    print(event)
    assert str(context.create_annotations) == event, f"Expected event: {event}, but got: {str(context.create_annotations)}"