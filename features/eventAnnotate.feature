Feature: Annotate events in splicing data
	To start with:
	- I have a splicing event with chrom, start, end, strand, and type
	- I have a transcript annotation file

	Firstly:
	- Annotate the starts
	- Annotate the ends
	- Find the intron overlapping the splicing event
	1. Does the start match? (canonical, cryptics, skipping)
	If yes: Does the end match? -> canonical
	If no: Send to end_matcher
	2. Does the end match? (cryptics, skipping)
	If yes: Send to start_matcher
	3. Is the start inside the intron?
	If yes: Send to end_matcher
	4. Is the end inside the intron?
	If yes: Send to start_matcher
	5. If all no: Send to an unmatched list

	Start Matcher:
	Matcher with starts

	End Matcher:
	Matcher with ends

	Matcher (coordinates, annotations, start/end):
	1. Canonical Splice-site
	2. Same Transcript Splice-site (skipping, multi-skipping)
	3. Different Transcript Splice-site (alternate-splicing, alternate-skipping)
	4. Unannotated Splice-site (cryptic, intronic)

	@wip
	Scenario Outline: Annotate events in cortar data using RefSeq MANE Select Annotations
		Given a splicing event with chrom <chrom>, start <start>, end <end>, strand <strand>, and type <type>
		And the transcript annotation is refseq_mane
		When the splicing events are annotated
		Then the result should be event <event> for transcript <transcript>

		Examples:
			| chrom | start    | end      | strand | transcript     | type | event                                               | intron | comment                      | comment2                                                                           |
			| chr9  | 34646787 | 34647088 | +      | NM_000155.4    | sj   | canonical exon 1-2 splicing                         | 1      | galt_intron_1_merged_sjir_4  | canonical +                                                                        |
			| chr9  | 34646787 | 34647088 | +      | NM_000155.4    | ir   | intron 1 retention                                  | 1      | galt_intron_1_merged_sjir_2  | intron retention +                                                                 |
			| chr9  | 34647259 | 34647831 | +      | NM_000155.4    | sj   | exon 3-4 skipping                                   | 2      | galt_intron_1_merged_sjir_9  | double exon skipping +                                                             |
			| chr9  | 34647259 | 34648114 | +      | NM_000155.4    | sj   | exon 3-4-5 skipping                                 | 2      | galt_intron_1_merged_sjir_6  | triple exon skipping +                                                             |
			| chr9  | 34646920 | 34647088 | +      | NM_000155.4    | sj   | intronic cryptic donor @ -133                       | 1      | galt_intron_1_merged_sjir_0  | cryptic donor exonic +                                                             |
			| chr9  | 34646787 | 34647120 | +      | NM_000155.4    | sj   | exonic cryptic acceptor @ +32                       | 1      | galt_intron_1_merged_sjir_5  | cryptic acceptor exonic +                                                          |
			| chr9  | 34647259 | 34647474 | +      | NM_000155.4    | sj   | intronic cryptic acceptor @ -17                     | 2      | galt_intron_1_merged_sjir_12 | cryptic acceptor intronic +                                                        |
			| chr9  | 34647259 | 34648091 | +      | NM_000155.4    | sj   | exon 3-4-5 skipping/intronic cryptic acceptor @ -23 | 2      | galt_intron_1_merged_sjir_7  | cryptic acceptor with skipping +                                                   |
			| chr9  | 34647259 | 34648047 | +      | NM_000155.4    | sj   | exon 3-4-5 skipping/intronic cryptic acceptor @ -67 | 2      | galt_intron_1_merged_sjir_8  | cryptic acceptor with skipping +                                                   |
			| chrX  | 40054043 | 40054255 | -      | NM_001123385.2 | sj   | canonical exon 13-14 splicing                       | 13     | galt_chrX_subset_BCOR        | canonical -                                                                        |
			| chrX  | 40054043 | 40054255 | -      | NM_001123385.2 | ir   | intron 13 retention                                 | 13     | galt_chrX_subset_BCOR        | intron retention -                                                                 |
			| chrX  | 40057322 | 40062745 | -      | NM_001123385.2 | sj   | exon 10 skipping                                    | 9      | galt_chrX_subset_BCOR        | single exon skipping -                                                             |
			| chrX  | 40077970 | 40177006 | -      | NM_001123385.2 | sj   | alternate donor exon 1-2 splicing (NM_001123383.1)  | 1      | galt_chrX_subset_BCOR        | alternate donor -                                                                  |
			| chrX  | 40081196 | 40081258 | +      | NM_001123385.2 | sj   | unannotated junctions                               | 1      | galt_chrX_subset_BCOR        | unannotated junctions + within intron -                                            |
			| chrX  | 40096604 | 40097680 | +      | NM_001123385.2 | sj   | unannotated junctions                               | 1      | galt_chrX_subset_BCOR        | unannotated junctions + within intron -                                            |
			| chrX  | 40097722 | 40098095 | +      | NM_001123385.2 | sj   | unannotated junctions                               | NA     | galt_chrX_subset_BCOR        | unannotated junctions + overlaps exon 1 -                                          |
			| chrX  | 40113978 | 40208244 | -      | NM_001123385.2 | sj   | unannotated junctions (intergenic)                  | NA     | galt_chrX_subset_BCOR        | unannotated junctions - within alternate intron 1 - & intergenic (mis-map warning) |
			| chrX  | 40115468 | 40117489 | +      | NM_001123385.2 | sj   | unannotated junctions (intergenic)                  | NA     | galt_chrX_subset_BCOR        | unannotated junctions + within alternate intron 1 - (genuine splicing?)            |
# |       |          |          | +      |                |                                                 |        |                              | alternate intronic acceptor + |
# |       |          |          | +      |                |                                                 |        |                              | alternate exonic donor +      |
# |       |          |          | +      |                |                                                 |        |                              | alternate exonic acceptor +   |
# |       |          |          | _      |                |                                                 |        |                              | double exon skipping -        |
# |       |          |          | _      |                |                                                 |        |                              | triple exon skipping -        |
# |       |          |          | _      |                |                                                 |        |                              | cryptic exonic donor -        |
# |       |          |          | _      |                |                                                 |        |                              | cryptic exonic acceptor -     |
# |       |          |          | _      |                |                                                 |        |                              | cryptic intronic donor -      |
# |       |          |          | _      |                |                                                 |        |                              | cryptic intronic acceptor -   |
# |       |          |          | _      |                |                                                 |        |                              | cryptic exonic acceptor -     |
# |       |          |          | _      |                |                                                 |        |                              | cryptic intronic acceptor -   |
# |       |          |          | _      |                |                                                 |        |                              | alternate intronic donor -    |
# |       |          |          | _      |                |                                                 |        |                              | alternate intronic acceptor - |
# |       |          |          | _      |                |                                                 |        |                              | alternate exonic donor -      |
# |       |          |          | _      |                |                                                 |        |                              | alternate exonic acceptor -   |
# |       |          |          | _      |                |                                                 |        |                              | intergenic -                  |
# |       |          |          | _      |                |                                                 |        |                              | gene fusion -                 |

