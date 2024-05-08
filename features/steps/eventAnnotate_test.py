from behave import given, when, then
from src.eventAnnotateProcessor import EventAnnotate



# Fetch the MANE transcript for a splicing event

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
    )
    context.transcript = context.annotation.get_mane_transcript()


@then(u'the result should be transcript {transcript}, gene {gene}, and warning {warning}')
def step_impl(context, transcript, gene, warning):
    print(context.transcript['transcript'])
    print(transcript)
    print(context.transcript['gene'])
    print(gene)
    print(context.transcript['warning'])
    print(warning)
    assert context.transcript['transcript'] == transcript
    assert context.transcript['gene'] == gene
    assert context.transcript['warning'] == warning
    


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
    )
    context.annotation.get_annotations(context.input['annotation'])
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    

@then(u'the resulting annotations of canonical splicing should be event {event} of type {event_type} at intron {intron}')
def step_impl(context, event, event_type, intron):
    print(str(context.create_annotations['event']))
    print(event)
    print(str(context.create_annotations['event_type']))
    print(event_type)
    print(str(context.create_annotations['introns']))
    print(intron)
    assert str(context.create_annotations['event']) == event, f"Expected event: {event}, but got: {str(context.create_annotations['event'])}"
    assert str(context.create_annotations['event_type']) == event_type, f"Expected event: {event_type}, but got: {str(context.create_annotations['event_type'])}"
    assert str(context.create_annotations['introns']) == intron, f"Expected event: {intron}, but got: {str(context.create_annotations['introns'])}"



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
    )
    context.annotation.get_annotations(context.input['annotation'])
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of exon skipping should be event {event} of type {event_type} at intron {intron}')
def step_impl(context, event, event_type, intron):
    print(str(context.create_annotations['event']))
    print(event)
    print(str(context.create_annotations['event_type']))
    print(event_type)
    print(str(context.create_annotations['introns']))
    print(intron)
    assert str(context.create_annotations['event']) == event, f"Expected event: {event}, but got: {str(context.create_annotations['event'])}"
    assert str(context.create_annotations['event_type']) == event_type, f"Expected event: {event_type}, but got: {str(context.create_annotations['event_type'])}"
    assert str(context.create_annotations['introns']) == intron, f"Expected event: {intron}, but got: {context.create_annotations['introns']}"



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
    )
    context.annotation.get_annotations(context.input['annotation'])
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of intron retention should be event {event} of type {event_type} at intron {intron}')
def step_impl(context, event, event_type, intron):
    print(str(context.create_annotations['event']))
    print(event)
    print(str(context.create_annotations['event_type']))
    print(event_type)
    print(str(context.create_annotations['introns']))
    print(intron)
    assert str(context.create_annotations['event']) == event, f"Expected event: {event}, but got: {str(context.create_annotations['event'])}"
    assert str(context.create_annotations['event_type']) == event_type, f"Expected event: {event_type}, but got: {str(context.create_annotations['event_type'])}"
    assert str(context.create_annotations['introns']) == intron, f"Expected event: {intron}, but got: {context.create_annotations['introns']}"



# CRYPTIC DONORS

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
    )
    context.annotation.get_annotations(context.input['annotation'])
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of cryptic donors should be event {event} of location {location} and type {event_type} at intron {intron} at a distance of {distance_from_authentic}')
def step_impl(context, event, event_type, intron, location, distance_from_authentic):
    print(str(context.create_annotations['event']))
    print(event)
    print(str(context.create_annotations['event_type']))
    print(event_type)
    print(str(context.create_annotations['introns']))
    print(intron)
    print(str(context.create_annotations['location']))
    print(location)
    print(str(context.create_annotations['distance_from_authentic']))
    print(distance_from_authentic)
    assert str(context.create_annotations['event']) == event, f"Expected event: {event}, but got: {str(context.create_annotations['event'])}"
    assert str(context.create_annotations['event_type']) == event_type, f"Expected event: {event_type}, but got: {str(context.create_annotations['event_type'])}"
    assert str(context.create_annotations['introns']) == intron, f"Expected event: {intron}, but got: {context.create_annotations['introns']}"
    assert str(context.create_annotations['location']) == location, f"Expected event: {location}, but got: {context.create_annotations['location']}"
    assert str(context.create_annotations['distance_from_authentic']) == distance_from_authentic, f"Expected event: {distance_from_authentic}, but got: {context.create_annotations['distance_from_authentic']}"



# CRYPTIC ACCEPTORS

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
    )
    context.annotation.get_annotations(context.input['annotation'])
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of cryptic acceptors should be event {event} of location {location} and type {event_type} at intron {intron} at a distance of {distance_from_authentic}')
def step_impl(context, event, event_type, intron, location, distance_from_authentic):
    print(str(context.create_annotations['event']))
    print(event)
    print(str(context.create_annotations['event_type']))
    print(event_type)
    print(str(context.create_annotations['introns']))
    print(intron)
    print(str(context.create_annotations['location']))
    print(location)
    print(str(context.create_annotations['distance_from_authentic']))
    print(distance_from_authentic)
    assert str(context.create_annotations['event']) == event, f"Expected event: {event}, but got: {str(context.create_annotations['event'])}"
    assert str(context.create_annotations['event_type']) == event_type, f"Expected event: {event_type}, but got: {str(context.create_annotations['event_type'])}"
    assert str(context.create_annotations['introns']) == intron, f"Expected event: {intron}, but got: {context.create_annotations['introns']}"
    assert str(context.create_annotations['location']) == location, f"Expected event: {location}, but got: {context.create_annotations['location']}"
    assert str(context.create_annotations['distance_from_authentic']) == distance_from_authentic, f"Expected event: {distance_from_authentic}, but got: {context.create_annotations['distance_from_authentic']}"



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
    )
    context.annotation.get_annotations(context.input['annotation'])
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of skip-cryps should be event {event} of location {location} and type {event_type} at intron {intron} at a distance of {distance_from_authentic}')
def step_impl(context, event, event_type, intron, location, distance_from_authentic):
    print(str(context.create_annotations['event']))
    print(event)
    print(str(context.create_annotations['event_type']))
    print(event_type)
    print(str(context.create_annotations['introns']))
    print(intron)
    print(str(context.create_annotations['location']))
    print(location)
    print(str(context.create_annotations['distance_from_authentic']))
    print(distance_from_authentic)
    assert str(context.create_annotations['event']) == event, f"Expected event: {event}, but got: {str(context.create_annotations['event'])}"
    assert str(context.create_annotations['event_type']) == event_type, f"Expected event: {event_type}, but got: {str(context.create_annotations['event_type'])}"
    assert str(context.create_annotations['introns']) == intron, f"Expected event: {intron}, but got: {context.create_annotations['introns']}"
    assert str(context.create_annotations['location']) == location, f"Expected event: {location}, but got: {context.create_annotations['location']}"
    assert str(context.create_annotations['distance_from_authentic']) == distance_from_authentic, f"Expected event: {distance_from_authentic}, but got: {context.create_annotations['distance_from_authentic']}"



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
    )
    context.annotation.get_annotations(context.input['annotation'])
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of unannotated junctions should be event {event} of location {location} and type {event_type} at intron {intron}')
def step_impl(context, event, event_type, intron, location):
    print(str(context.create_annotations['event']))
    print(event)
    print(str(context.create_annotations['event_type']))
    print(event_type)
    print(str(context.create_annotations['introns']))
    print(intron)
    assert str(context.create_annotations['event']) == event, f"Expected event: {event}, but got: {str(context.create_annotations['event'])}"
    assert str(context.create_annotations['event_type']) == event_type, f"Expected event: {event_type}, but got: {str(context.create_annotations['event_type'])}"
    assert str(context.create_annotations['introns']) == intron, f"Expected event: {intron}, but got: {context.create_annotations['introns']}"



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
    )
    context.annotation.get_annotations(context.input['annotation'])
    context.start = context.annotation.reference_match('start')
    context.end = context.annotation.reference_match('end')
    
    print(context.input["start"])
    print(context.input["end"])
    print(context.start)
    print(context.end)
    
    context.create_annotations = context.annotation.fetch_transcript_annotations(context.start, context.end)
    
@then(u'the resulting annotations of alternate splicing should be event {event} of location {location} and type {event_type} at intron {intron} at a distance of {distance_from_authentic}')
def step_impl(context, event, event_type, intron, location, distance_from_authentic):
    print(str(context.create_annotations['event']))
    print(event)
    print(str(context.create_annotations['event_type']))
    print(event_type)
    print(str(context.create_annotations['introns']))
    print(intron)
    print(str(context.create_annotations['location']))
    print(location)
    print(str(context.create_annotations['distance_from_authentic']))
    print(distance_from_authentic)
    assert str(context.create_annotations['event']) == event, f"Expected event: {event}, but got: {str(context.create_annotations['event'])}"
    assert str(context.create_annotations['event_type']) == event_type, f"Expected event: {event_type}, but got: {str(context.create_annotations['event_type'])}"
    assert str(context.create_annotations['introns']) == intron, f"Expected event: {intron}, but got: {context.create_annotations['introns']}"
    assert str(context.create_annotations['location']) == location, f"Expected event: {location}, but got: {context.create_annotations['location']}"
    assert str(context.create_annotations['distance_from_authentic']) == distance_from_authentic, f"Expected event: {distance_from_authentic}, but got: {context.create_annotations['distance_from_authentic']}"
