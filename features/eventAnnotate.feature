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
			| test | chrom | start    | end      | strand | type | transcript     | event                                                  | intron | comment                      | comment2                                                                                                         |
			| 1    | chr9  | 34646787 | 34647088 | +      | sj   | NM_000155.4    | canonical exon 1-2 splicing                            | 1      | galt_intron_1_merged_sjir_4  | canonical +                                                                                                      |
			| 2    | chr9  | 34646787 | 34647088 | +      | ir   | NM_000155.4    | intron 1 retention                                     | 1      | galt_intron_1_merged_sjir_2  | intron retention +                                                                                               |
			| 3    | chr9  | 34647259 | 34647831 | +      | sj   | NM_000155.4    | exon 3-4 skipping                                      | 2      | galt_intron_1_merged_sjir_9  | double exon skipping +                                                                                           |
			| 4    | chr9  | 34647259 | 34648114 | +      | sj   | NM_000155.4    | exon 3-4-5 skipping                                    | 2      | galt_intron_1_merged_sjir_6  | triple exon skipping +                                                                                           |
			| 5    | chr9  | 34646920 | 34647088 | +      | sj   | NM_000155.4    | intronic cryptic donor @ +133                          | 1      | galt_intron_1_merged_sjir_0  | cryptic donor exonic +                                                                                           |
			| 6    | chr9  | 34646787 | 34647120 | +      | sj   | NM_000155.4    | exonic cryptic acceptor @ +32                          | 1      | galt_intron_1_merged_sjir_5  | cryptic acceptor exonic +                                                                                        |
			| 7    | chr9  | 34647259 | 34647474 | +      | sj   | NM_000155.4    | intronic cryptic acceptor @ -18                        | 2      | galt_intron_1_merged_sjir_12 | cryptic acceptor intronic +                                                                                      |
			| 8    | chr9  | 34647259 | 34648091 | +      | sj   | NM_000155.4    | exon 3-4-5 skipping/intronic cryptic acceptor @ -24    | 2      | galt_intron_1_merged_sjir_7  | cryptic acceptor with skipping +                                                                                 |
			| 9    | chr9  | 34647259 | 34648047 | +      | sj   | NM_000155.4    | exon 3-4-5 skipping/intronic cryptic acceptor @ -68    | 2      | galt_intron_1_merged_sjir_8  | cryptic acceptor with skipping +                                                                                 |
			| 10   | chrX  | 40054043 | 40054255 | -      | sj   | NM_001123385.2 | canonical exon 13-14 splicing                          | 13     | galt_chrX_subset_BCOR        | canonical -                                                                                                      |
			| 11   | chrX  | 40054043 | 40054255 | -      | ir   | NM_001123385.2 | intron 13 retention                                    | 13     | galt_chrX_subset_BCOR        | intron retention -                                                                                               |
			| 12   | chrX  | 40057322 | 40062745 | -      | sj   | NM_001123385.2 | exon 10 skipping                                       | 9      | galt_chrX_subset_BCOR        | single exon skipping -                                                                                           |
			| 13   | chrX  | 40081196 | 40081258 | +      | sj   | NM_001123385.2 | unannotated intronic junction (opposite strand)        | 1      | galt_chrX_subset_BCOR        | unannotated junctions + (NM_001123385.2) within intron -                                                         |
			| 14   | chrX  | 40096604 | 40097680 | +      | sj   | NM_001123385.2 | unannotated intronic/exonic junction (opposite strand) | 1      | galt_chrX_subset_BCOR        | unannotated junctions + (NM_001123385.2) within intron -                                                         |
			| 15   | chrX  | 40097722 | 40098095 | +      | sj   | NM_001123385.2 | unknown intergenic/exonic junction (opposite strand)   | NA     | galt_chrX_subset_BCOR        | unannotated junctions + (NM_001123385.2) overlaps exon 1 -                                                       |
			| 16   | chrX  | 40113978 | 40208244 | -      | sj   | unknown        | unknown junction                                       | NA     | galt_chrX_subset_BCOR        | unannotated junctions - (NM_001123385.2) within alternate intron 1 - & intergenic (mis-map warning)              |
			| 17   | chrX  | 40115468 | 40117489 | +      | sj   | unknown        | unknown junction                                       | NA     | galt_chrX_subset_BCOR        | unannotated junctions + (NM_001123385.2) within alternate intron 1 - (genuine splicing?)                         |
			| 18   | chrX  | 40077970 | 40177006 | -      | sj   | NM_001123385.2 | intergenic cryptic donor @ -79792                      | 1      | galt_chrX_subset_BCOR        | alternate donor exon 1-2 splicing (NM_001123383.1) alternate donor -                                             |
			| 19   | chrX  | 41086110 | 41123470 | +      | sj   | NM_001039591.3 | canonical exon 1-2 splicing                            | 1      | galt_chrX_subset_USP9X       | canonical +                                                                                                      |
			| 20   | chrX  | 41087666 | 41123470 | +      | sj   | NM_001039591.3 | intronic cryptic donor @ +1556                         | 1      | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 21   | chrX  | 41101571 | 41102912 | -      | sj   | NM_001039591.3 | unannotated intronic junction (opposite strand)        | 1      | galt_chrX_subset_USP9X       | 2                                                                                                                |
			| 22   | chrX  | 41114762 | 41123470 | +      | sj   | NM_001039591.3 | intronic cryptic donor @ +28652                        | 1      | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 23   | chrX  | 41123725 | 41128999 | +      | sj   | NM_001039591.3 | canonical exon 2-3 splicing                            | 2      | galt_chrX_subset_USP9X       | 2                                                                                                                |
			| 24   | chrX  | 41134838 | 41140655 | +      | sj   | NM_001039591.3 | exon 6 skipping                                        |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 25   | chrX  | 41141436 | 41143290 | +      | sj   | NM_001039591.3 | intronic cryptic donor @ +4                            |        | galt_chrX_subset_USP9X       | super interesting, uses an AT acceptor at -4 but soft clipped to +4 donor position, not a canonical donor site 1 |
			| 26   | chrX  | 41144347 | 41144521 | +      | sj   | NM_001039591.3 | intronic cryptic donor @ +903                          |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 27   | chrX  | 41144627 | 41146305 | +      | sj   | NM_001039591.3 | intronic cryptic acceptor @ -2064                      |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 28   | chrX  | 41144627 | 41148373 | +      | sj   | NM_001039591.3 | exonic cryptic acceptor @ +5                           |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 29   | chrX  | 41148576 | 41152932 | +      | sj   | NM_001039591.3 | exon 13 skipping/intronic cryptic acceptor @ -16       |        | galt_chrX_subset_USP9X       | cryptic annotated in other transcripts but not with skipping                                                     |
			| 30   | chrX  | 41149755 | 41151932 | +      | sj   | NM_001039591.3 | unannotated intronic junction                          |        | galt_chrX_subset_USP9X       | this one also skips exon 13                                                                                      |
			| 31   | chrX  | 41153082 | 41167481 | +      | sj   | NM_001039591.3 | exon 15-16 skipping                                    |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 32   | chrX  | 41198618 | 41205302 | +      | sj   | NM_001039591.3 | exon 31 skipping/exonic cryptic donor @ -133           |        | galt_chrX_subset_USP9X       | actually an exonic cryptic donor plus skipping of exon 31                                                        |
			| 33   | chrX  | 41216059 | 41217219 | +      | sj   | NM_001039591.3 | exonic cryptic donor @ +5914   (wrong)                 |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 34   | chrX  | 41223403 | 41224806 | +      | sj   | NM_001039591.3 | exonic cryptic acceptor @ +69 (wrong)                  |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 35   | chrX  | 41229636 | 41230503 | +      | sj   | NM_001039591.3 | unannotated exonic junctions (wrong)                   |        | galt_chrX_subset_USP9X       | 0                                                                                                                |



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

