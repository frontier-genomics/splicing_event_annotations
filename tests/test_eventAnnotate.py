import pytest
from pytest_bdd import given, when, then, scenarios, parsers
from splicing_event_annotator.eventAnnotateProcessor import EventAnnotate
import logging

logging.basicConfig(level=logging.WARNING)

scenarios('eventAnnotate.feature')

# Fetch the MANE transcript for a splicing event

@given(parsers.parse('a splicing event with chrom {chrom}, start {start}, end {end}, and strand {strand}'), target_fixture='event_input')
def splicing_event(chrom, start, end, strand):
    """Create splicing event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': "NA",
        'type': ""
    }

@when('the overlapping transcripts are fetched', target_fixture='transcript_result')
def fetch_transcripts(event_input, refgene_annotations):
    """Fetch overlapping transcripts."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation.get_mane_transcript()

@then(parsers.parse('the result should be transcript {transcript}, gene {gene}, and warning {warning}'))
def verify_transcript(transcript_result, transcript, gene, warning):
    """Verify transcript result."""
    print(transcript_result['transcript'])
    print(transcript)
    print(transcript_result['gene'])
    print(gene)
    print(transcript_result['warning'])
    print(warning)
    assert transcript_result['transcript'] == transcript
    assert transcript_result['gene'] == gene
    assert transcript_result['warning'] == warning

# CANONICAL SPLICING

@given(parsers.parse('a canonical splicing event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}'), target_fixture='event_input')
def canonical_splicing_event(chrom, start, end, strand, transcript, type):
    """Create canonical splicing event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': transcript,
        'type': type
    }

@given(parsers.parse('the annotation dataset to be used for annotating canonical splicing is {annotation}'))
def canonical_annotation_dataset(annotation):
    """Store annotation dataset for canonical splicing."""
    pass  # Dataset is stored but not used in this step

@when('the canonical splicing events are annotated with the annotation dataset', target_fixture='canonical_annotation')
def annotate_canonical(event_input, refgene_annotations):
    """Annotate canonical splicing events."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation

@then(parsers.parse('the resulting annotations of canonical splicing should be event {event} of type {event_type} at intron {intron}'))
def verify_canonical_annotation(canonical_annotation, event, event_type, intron):
    """Verify canonical splicing annotations."""
    print(str(canonical_annotation.event))
    print(event)
    print(str(canonical_annotation.event_type))
    print(event_type)
    print(str(canonical_annotation.introns))
    print(intron)
    assert str(canonical_annotation.event) == event, f"Expected event: {event}, but got: {str(canonical_annotation.event)}"
    assert str(canonical_annotation.event_type) == event_type, f"Expected event: {event_type}, but got: {str(canonical_annotation.event_type)}"
    assert str(canonical_annotation.introns) == intron, f"Expected event: {intron}, but got: {str(canonical_annotation.introns)}"

# EXON SKIPPING

@given(parsers.parse('an exon skipping event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}'), target_fixture='event_input')
def exon_skipping_event(chrom, start, end, strand, transcript, type):
    """Create exon skipping event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': transcript,
        'type': type
    }

@given(parsers.parse('the annotation dataset to be used for annotating exon skipping is {annotation}'))
def exon_skipping_annotation_dataset(annotation):
    """Store annotation dataset for exon skipping."""
    pass

@when('the exon skipping events are annotated with the annotation dataset', target_fixture='exon_skipping_annotation')
def annotate_exon_skipping(event_input, refgene_annotations):
    """Annotate exon skipping events."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation

@then(parsers.parse('the resulting annotations of exon skipping should be event {event} of type {event_type} at intron {intron}'))
def verify_exon_skipping_annotation(exon_skipping_annotation, event, event_type, intron):
    """Verify exon skipping annotations."""
    print(str(exon_skipping_annotation.event))
    print(event)
    print(str(exon_skipping_annotation.event_type))
    print(event_type)
    print(str(exon_skipping_annotation.introns))
    print(intron)
    assert str(exon_skipping_annotation.event) == event, f"Expected event: {event}, but got: {str(exon_skipping_annotation.event)}"
    assert str(exon_skipping_annotation.event_type) == event_type, f"Expected event: {event_type}, but got: {str(exon_skipping_annotation.event_type)}"
    assert str(exon_skipping_annotation.introns) == intron, f"Expected event: {intron}, but got: {exon_skipping_annotation.introns}"

# INTRON RETENTION

@given(parsers.parse('an intron retention event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}'), target_fixture='event_input')
def intron_retention_event(chrom, start, end, strand, transcript, type):
    """Create intron retention event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': transcript,
        'type': type
    }

@given(parsers.parse('the annotation dataset to be used for annotating intron retention is {annotation}'))
def intron_retention_annotation_dataset(annotation):
    """Store annotation dataset for intron retention."""
    pass

@when('the intron retention events are annotated with the annotation dataset', target_fixture='intron_retention_annotation')
def annotate_intron_retention(event_input, refgene_annotations):
    """Annotate intron retention events."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation

@then(parsers.parse('the resulting annotations of intron retention should be event {event} of type {event_type} at intron {intron}'))
def verify_intron_retention_annotation(intron_retention_annotation, event, event_type, intron):
    """Verify intron retention annotations."""
    print(str(intron_retention_annotation.event))
    print(event)
    print(str(intron_retention_annotation.event_type))
    print(event_type)
    print(str(intron_retention_annotation.introns))
    print(intron)
    assert str(intron_retention_annotation.event) == event, f"Expected event: {event}, but got: {str(intron_retention_annotation.event)}"
    assert str(intron_retention_annotation.event_type) == event_type, f"Expected event: {event_type}, but got: {str(intron_retention_annotation.event_type)}"
    assert str(intron_retention_annotation.introns) == intron, f"Expected event: {intron}, but got: {intron_retention_annotation.introns}"

# CRYPTIC DONORS

@given(parsers.parse('a cryptic donor event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}'), target_fixture='event_input')
def cryptic_donor_event(chrom, start, end, strand, transcript, type):
    """Create cryptic donor event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': transcript,
        'type': type
    }

@given(parsers.parse('the annotation dataset to be used for annotating cryptic donors is {annotation}'))
def cryptic_donor_annotation_dataset(annotation):
    """Store annotation dataset for cryptic donors."""
    pass

@when('the cryptic donor events are annotated with the annotation dataset', target_fixture='cryptic_donor_annotation')
def annotate_cryptic_donor(event_input, refgene_annotations):
    """Annotate cryptic donor events."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation

@then(parsers.parse('the resulting annotations of cryptic donors should be event {event} of location {location} and type {event_type} at intron {intron} at a distance of {distance_from_authentic}'))
def verify_cryptic_donor_annotation(cryptic_donor_annotation, event, event_type, intron, location, distance_from_authentic):
    """Verify cryptic donor annotations."""
    print(str(cryptic_donor_annotation.event))
    print(event)
    print(str(cryptic_donor_annotation.event_type))
    print(event_type)
    print(str(cryptic_donor_annotation.introns))
    print(intron)
    print(str(cryptic_donor_annotation.location))
    print(location)
    print(str(cryptic_donor_annotation.distance_from_authentic))
    print(distance_from_authentic)
    assert str(cryptic_donor_annotation.event) == event, f"Expected event: {event}, but got: {str(cryptic_donor_annotation.event)}"
    assert str(cryptic_donor_annotation.event_type) == event_type, f"Expected event: {event_type}, but got: {str(cryptic_donor_annotation.event_type)}"
    assert str(cryptic_donor_annotation.introns) == intron, f"Expected event: {intron}, but got: {cryptic_donor_annotation.introns}"
    assert str(cryptic_donor_annotation.location) == location, f"Expected event: {location}, but got: {cryptic_donor_annotation.location}"
    assert str(cryptic_donor_annotation.distance_from_authentic) == distance_from_authentic, f"Expected event: {distance_from_authentic}, but got: {cryptic_donor_annotation.distance_from_authentic}"

# CRYPTIC ACCEPTORS

@given(parsers.parse('a cryptic acceptor event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}'), target_fixture='event_input')
def cryptic_acceptor_event(chrom, start, end, strand, transcript, type):
    """Create cryptic acceptor event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': transcript,
        'type': type
    }

@given(parsers.parse('the annotation dataset to be used for annotating cryptic acceptors is {annotation}'))
def cryptic_acceptor_annotation_dataset(annotation):
    """Store annotation dataset for cryptic acceptors."""
    pass

@when('the cryptic acceptor events are annotated with the annotation dataset', target_fixture='cryptic_acceptor_annotation')
def annotate_cryptic_acceptor(event_input, refgene_annotations):
    """Annotate cryptic acceptor events."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation

@then(parsers.parse('the resulting annotations of cryptic acceptors should be event {event} of location {location} and type {event_type} at intron {intron} at a distance of {distance_from_authentic}'))
def verify_cryptic_acceptor_annotation(cryptic_acceptor_annotation, event, event_type, intron, location, distance_from_authentic):
    """Verify cryptic acceptor annotations."""
    print(str(cryptic_acceptor_annotation.event))
    print(event)
    print(str(cryptic_acceptor_annotation.event_type))
    print(event_type)
    print(str(cryptic_acceptor_annotation.introns))
    print(intron)
    print(str(cryptic_acceptor_annotation.location))
    print(location)
    print(str(cryptic_acceptor_annotation.distance_from_authentic))
    print(distance_from_authentic)
    assert str(cryptic_acceptor_annotation.event) == event, f"Expected event: {event}, but got: {str(cryptic_acceptor_annotation.event)}"
    assert str(cryptic_acceptor_annotation.event_type) == event_type, f"Expected event: {event_type}, but got: {str(cryptic_acceptor_annotation.event_type)}"
    assert str(cryptic_acceptor_annotation.introns) == intron, f"Expected event: {intron}, but got: {cryptic_acceptor_annotation.introns}"
    assert str(cryptic_acceptor_annotation.location) == location, f"Expected event: {location}, but got: {cryptic_acceptor_annotation.location}"
    assert str(cryptic_acceptor_annotation.distance_from_authentic) == distance_from_authentic, f"Expected event: {distance_from_authentic}, but got: {cryptic_acceptor_annotation.distance_from_authentic}"

# SKIPPING WITH A CRYPTIC (SKIP CRYPS)

@given(parsers.parse('a skip-cryp event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}'), target_fixture='event_input')
def skip_cryp_event(chrom, start, end, strand, transcript, type):
    """Create skip-cryp event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': transcript,
        'type': type
    }

@given(parsers.parse('the annotation dataset to be used for annotating skip-cryps is {annotation}'))
def skip_cryp_annotation_dataset(annotation):
    """Store annotation dataset for skip-cryps."""
    pass

@when('the skip-cryp events are annotated with the annotation dataset', target_fixture='skip_cryp_annotation')
def annotate_skip_cryp(event_input, refgene_annotations):
    """Annotate skip-cryp events."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation

@then(parsers.parse('the resulting annotations of skip-cryps should be event {event} of location {location} and type {event_type} at intron {intron} at a distance of {distance_from_authentic}'))
def verify_skip_cryp_annotation(skip_cryp_annotation, event, event_type, intron, location, distance_from_authentic):
    """Verify skip-cryp annotations."""
    print(str(skip_cryp_annotation.event))
    print(event)
    print(str(skip_cryp_annotation.event_type))
    print(event_type)
    print(str(skip_cryp_annotation.introns))
    print(intron)
    print(str(skip_cryp_annotation.location))
    print(location)
    print(str(skip_cryp_annotation.distance_from_authentic))
    print(distance_from_authentic)
    assert str(skip_cryp_annotation.event) == event, f"Expected event: {event}, but got: {str(skip_cryp_annotation.event)}"
    assert str(skip_cryp_annotation.event_type) == event_type, f"Expected event: {event_type}, but got: {str(skip_cryp_annotation.event_type)}"
    assert str(skip_cryp_annotation.introns) == intron, f"Expected event: {intron}, but got: {skip_cryp_annotation.introns}"
    assert str(skip_cryp_annotation.location) == location, f"Expected event: {location}, but got: {skip_cryp_annotation.location}"
    assert str(skip_cryp_annotation.distance_from_authentic) == distance_from_authentic, f"Expected event: {distance_from_authentic}, but got: {skip_cryp_annotation.distance_from_authentic}"

# UNANNOTATED JUNCTIONS

@given(parsers.parse('an unannotated junction with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}'), target_fixture='event_input')
def unannotated_junction_event(chrom, start, end, strand, transcript, type):
    """Create unannotated junction event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': transcript,
        'type': type
    }

@given(parsers.parse('the annotation dataset to be used for annotating unannotated junctions is {annotation}'))
def unannotated_junction_annotation_dataset(annotation):
    """Store annotation dataset for unannotated junctions."""
    pass

@when('the unannotated junctions are annotated with the annotation dataset', target_fixture='unannotated_junction_annotation')
def annotate_unannotated_junction(event_input, refgene_annotations):
    """Annotate unannotated junction events."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation

@then(parsers.parse('the resulting annotations of unannotated junctions should be event {event} of location {location} and type {event_type} at intron {intron}'))
def verify_unannotated_junction_annotation(unannotated_junction_annotation, event, event_type, intron, location):
    """Verify unannotated junction annotations."""
    print(str(unannotated_junction_annotation.event))
    print(event)
    print(str(unannotated_junction_annotation.event_type))
    print(event_type)
    print(str(unannotated_junction_annotation.introns))
    print(intron)
    assert str(unannotated_junction_annotation.event) == event, f"Expected event: {event}, but got: {str(unannotated_junction_annotation.event)}"
    assert str(unannotated_junction_annotation.event_type) == event_type, f"Expected event: {event_type}, but got: {str(unannotated_junction_annotation.event_type)}"
    assert str(unannotated_junction_annotation.introns) == intron, f"Expected event: {intron}, but got: {unannotated_junction_annotation.introns}"

# ALTERNATE SPLICING

@given(parsers.parse('an alternate splicing event with chrom {chrom}, start {start}, end {end}, strand {strand}, transcript {transcript}, and type {type}'), target_fixture='event_input')
def alternate_splicing_event(chrom, start, end, strand, transcript, type):
    """Create alternate splicing event input."""
    return {
        'chrom': chrom,
        'start': start,
        'end': end,
        'strand': strand,
        'transcript': transcript,
        'type': type
    }

@given(parsers.parse('the annotation dataset to be used for annotating alternate splicing is {annotation}'))
def alternate_splicing_annotation_dataset(annotation):
    """Store annotation dataset for alternate splicing."""
    pass

@when('the alternate splicing events are annotated with the annotation dataset', target_fixture='alternate_splicing_annotation')
def annotate_alternate_splicing(event_input, refgene_annotations):
    """Annotate alternate splicing events."""
    annotation = EventAnnotate(
        chrom=event_input['chrom'],
        start=event_input['start'],
        end=event_input['end'],
        strand=event_input['strand'],
        transcript=event_input['transcript'],
        input_type=event_input['type']
    )
    annotation.process('refseq', 'hg38', get_annotations=refgene_annotations)
    return annotation

@then(parsers.parse('the resulting annotations of alternate splicing should be event {event} of location {location} and type {event_type} at intron {intron} at a distance of {distance_from_authentic}'))
def verify_alternate_splicing_annotation(alternate_splicing_annotation, event, event_type, intron, location, distance_from_authentic):
    """Verify alternate splicing annotations."""
    print(str(alternate_splicing_annotation.event))
    print(event)
    print(str(alternate_splicing_annotation.event_type))
    print(event_type)
    print(str(alternate_splicing_annotation.introns))
    print(intron)
    print(str(alternate_splicing_annotation.location))
    print(location)
    print(str(alternate_splicing_annotation.distance_from_authentic))
    print(distance_from_authentic)
    assert str(alternate_splicing_annotation.event) == event, f"Expected event: {event}, but got: {str(alternate_splicing_annotation.event)}"
    assert str(alternate_splicing_annotation.event_type) == event_type, f"Expected event: {event_type}, but got: {str(alternate_splicing_annotation.event_type)}"
    assert str(alternate_splicing_annotation.introns) == intron, f"Expected event: {intron}, but got: {alternate_splicing_annotation.introns}"
    assert str(alternate_splicing_annotation.location) == location, f"Expected event: {location}, but got: {alternate_splicing_annotation.location}"
    assert str(alternate_splicing_annotation.distance_from_authentic) == distance_from_authentic, f"Expected event: {distance_from_authentic}, but got: {alternate_splicing_annotation.distance_from_authentic}"
