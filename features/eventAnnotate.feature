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
			| chrom | start    | end      | strand | transcript     | event                                           | intron | comment                      | comment2                    |
			| chr9  | 34646787 | 34647088 | +      | NM_000155.4    | canonical 1-2 splicing                          | 1      | galt_intron_1_merged_sjir_4  | canonical +                 |
			| chr9  | 34646787 | 34647088 | +      | NM_000155.4    | intron 1 retention                              | 1      | galt_intron_1_merged_sjir_2  | intron retention +          |
			| chr9  | 34647259 | 34647831 | +      | NM_000155.4    | exon 3-4 skipping                               | 2      | galt_intron_1_merged_sjir_9  | double exon skipping +      |
			| chr9  | 34647259 | 34648114 | +      | NM_000155.4    | exon 3-4-5 skipping                             | 2      | galt_intron_1_merged_sjir_6  | triple exon skipping +      |
			| chr9  | 34646920 | 34647088 | +      | NM_000155.4    | exonic cryptic donor @ -133                     | 1      | galt_intron_1_merged_sjir_0  | cryptic donor exonic +      |
			| chr9  | 34646787 | 34647120 | +      | NM_000155.4    | exonic cryptic acceptor @ +32                   | 1      | galt_intron_1_merged_sjir_5  | cryptic acceptor exonic +   |
			| chr9  | 34647259 | 34647474 | +      | NM_000155.4    | intronic cryptic acceptor @ -17                 | 2      | galt_intron_1_merged_sjir_12 | cryptic acceptor intronic + |
			| chr9  | 34647259 | 34648091 | +      | NM_000155.4    | cryptic acceptor                                | 2      | galt_intron_1_merged_sjir_7  | funky 1                     |
			| chr9  | 34647259 | 34648047 | +      | NM_000155.4    | cryptic acceptor                                | 2      | galt_intron_1_merged_sjir_8  | funky 2                     |
# |       |          |          | +      |                |                                                 |        |                              | alternate intronic acceptor + |
# |       |          |          | +      |                |                                                 |        |                              | alternate exonic donor +      |
# |       |          |          | +      |                |                                                 |        |                              | alternate exonic acceptor +   |
# |       |          |          | _      |                |                                                 |        |                              | canonical -                   |
# |       |          |          | _      |                |                                                 |        |                              | intron retention -            |
# |       |          |          | _      |                |                                                 |        |                              | single exon skipping -        |
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

