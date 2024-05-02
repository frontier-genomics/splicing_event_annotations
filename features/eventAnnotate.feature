Feature: Annotate events in cortar data using RefSeq Curated Annotations

	Scenario Outline: Canonical Splicing
		Given a canonical splicing event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating canonical splicing is refseq_curated
		When the canonical splicing events are annotated with the annotation dataset
		Then the resulting annotations of canonical splicing should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event                         | intron | comment |
			| chr9  | 34646787 | 34647088 | +      | sj   | NM_000155.4    | canonical exon 1-2 splicing   | 1      |         |
			| chrX  | 40054043 | 40054255 | -      | sj   | NM_001123385.2 | canonical exon 13-14 splicing | 13     |         |
			| chrX  | 41086110 | 41123470 | +      | sj   | NM_001039591.3 | canonical exon 1-2 splicing   | 1      |         |
			| chrX  | 41123725 | 41128999 | +      | sj   | NM_001039591.3 | canonical exon 2-3 splicing   | 2      |         |
			| chrX  | 13749534 | 13751248 | +      | sj   | NM_003611.3    | canonical exon 9-10 splicing  | 9      |         |
			| chrX  | 13751369 | 13753367 | +      | sj   | NM_003611.3    | canonical exon 10-11 splicing | 10     |         |
			| chrX  | 13749534 | 13751248 | +      | sj   | NM_001330210.2 | canonical exon 10-11 splicing | 10     |         |
			| chrX  | 13751369 | 13753367 | +      | sj   | NM_001330210.2 | canonical exon 11-12 splicing | 11     |         |
			| chrX  | 13749534 | 13753367 | +      | sj   | NM_001330209.2 | canonical exon 9-10 splicing  | 9      |         |
			| chrX  | 53256062 | 53291894 | -      | sj   | NM_001111125.3 | canonical exon 2-3 splicing   | 2      |         |
			| chrX  | 53256062 | 53291894 | -      | sj   | NM_001410736.1 | canonical exon 2-3 splicing   | 2      |         |
			| chrX  | 53279660 | 53281495 | -      | sj   | NM_001243197.2 | canonical exon 1-2 splicing   | 1      |         |
			| chrX  | 53267065 | 53279563 | -      | sj   | NM_001243197.2 | canonical exon 2-3 splicing   | 2      |         |
			| chrX  | 53279660 | 53281495 | -      | sj   | NM_015075.2    | canonical exon 1-2 splicing   | 1      |         |
			| chrX  | 53256062 | 53279563 | -      | sj   | NM_015075.2    | canonical exon 2-3 splicing   | 2      |         |


	Scenario Outline: Alternate Splicing
		Given an alternate splicing event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating alternate splicing is refseq_curated
		When the alternate splicing events are annotated with the annotation dataset
		Then the resulting annotations of alternate splicing should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event                                                                | intron | comment                                            |
			| chrX  | 13749534 | 13753367 | +      | sj   | NM_003611.3    | alternate exon 10 skipping (NM_001330209.2 exon 9-10)                |        | single exon skipping; 1 alt transcript             |
			| chrX  | 13749534 | 13753367 | +      | sj   | NM_001330210.2 | alternate exon 11 skipping (NM_001330209.2 exon 9-10)                |        | single exon skipping; 1 alt transcript             |
			| chr9  | 34647259 | 34647831 | +      | sj   | NM_000155.4    | alternate exon 3-4 skipping (NM_001258332.2 exon 2-3)                | 2      | double exon skipping; 1 alt transcript             |
			| chrX  | 41229828 | 41230500 | +      | sj   | NM_001039591.3 | alternate intron 43 donor @ +49 (NM_001039590.3, NM_001410748.1)     | 43     | intronic donor; 2 alt transcripts                  |
			| chrX  | 41229828 | 41230500 | +      | sj   | NM_001410749.1 | alternate intron 44 donor @ +49 (NM_001039590.3, NM_001410748.1)     | 44     | intronic donor; 2 alt transcripts                  |
			| chrX  | 13751369 | 13753367 | +      | sj   | NM_001330209.2 | alternate intron 9 donor @ +1836 (NM_001330210.2, NM_003611.3)       | 9      | intronic donor; 2 alt transcripts                  |
			| chrX  | 53256062 | 53279563 | -      | sj   | NM_001111125.3 | alternate intron 2 donor @ +12331 (NM_015075.2 exon 2-3)             | 2      | intronic donor; 1 alt transcript                   |
			| chrX  | 53256062 | 53279563 | -      | sj   | NM_001410736.1 | alternate intron 2 donor @ +12331 (NM_015075.2 exon 2-3)             |        | intronic donor; 1 alt transcript                   |
			| chrX  | 41087666 | 41123470 | +      | sj   | NM_001039591.3 | alternate intron 1 donor @ +1557 (NM_001410748.1, NM_001410749.1)    | 1      | intronic donor; 2 alt transcripts                  |
			| chrX  | 13749534 | 13751248 | +      | sj   | NM_001330209.2 | alternate intron 9 acceptor @ -2120 (NM_001330210.2, NM_003611.3)    | 9      | intronic acceptor; 2 alt transcripts               |
			| chrX  | 53267065 | 53279563 | -      | sj   | NM_015075.2    | alternate intron 2 acceptor @ -11004 (NM_001243197.2 exon 2-3)       |        | intronic acceptor; 1 alt transcript                |
			| chrX  | 53256062 | 53291894 | -      | sj   | NM_015075.2    | alternate intergenic donor @ -12331 (NM_001111125.3, NM_001410736.1) |        | intergenic donor; 2 alt transcripts                |
			| chrX  | 40077970 | 40177006 | -      | sj   | NM_001123385.2 | alternate intergenic donor @ -79792 (NM_001123383.1, NM_001123384.2) | 1      | intergenic donor; 2 alt transcripts                |
			| chrX  | 53256062 | 53279563 | -      | sj   | NM_001243197.2 | alternate intergenic acceptor @ +11003 (NM_015075.2 exon 2-3)        |        | intergenic acceptor; 1 alt transcript              |
			| chrX  | 53279660 | 53281495 | -      | sj   | NM_001111125.3 | alternate intronic junction (NM_001243197.2, NM_015075.2)            |        | unannotated intronic junction; 2 alt transcripts   |
			| chrX  | 53267065 | 53279563 | -      | sj   | NM_001111125.3 | alternate intronic junction (NM_001243197.2 exon 2-3)                |        | unannotated intronic junction; 1 alt transcript    |
			| chrX  | 53279660 | 53281495 | -      | sj   | NM_001410736.1 | alternate intronic junction (NM_001243197.2, NM_015075.2)            |        | unannotated intronic junction; 2 alt transcripts   |
			| chrX  | 53267065 | 53279563 | -      | sj   | NM_001410736.1 | alternate intronic junction (NM_001243197.2 exon 2-3)                |        | unannotated intronic junction; 1 alt transcript    |
			| chrX  | 53256062 | 53291894 | -      | sj   | NM_001243197.2 | alternate intergenic junction (NM_001111125.3, NM_001410736.1)       |        | unannotated intergenic junction; 2 alt transcripts |


	Scenario Outline: Exon Skipping
		Given an exon skipping event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating exon skipping is refseq_curated
		When the exon skipping events are annotated with the annotation dataset
		Then the resulting annotations of exon skipping should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event               | intron     | comment |
			| chrX  | 40057322 | 40062745 | -      | sj   | NM_001123385.2 | exon 10 skipping    | 9, 10      | single  |
			| chrX  | 41134838 | 41140655 | +      | sj   | NM_001039591.3 | exon 6 skipping     | 5, 6       | single  |
			| chrX  | 41153082 | 41167481 | +      | sj   | NM_001039591.3 | exon 15-16 skipping | 14, 15, 16 | double  |
			| chr9  | 34647259 | 34648114 | +      | sj   | NM_000155.4    | exon 3-4-5 skipping | 2, 3, 4, 5 | triple  |


	Scenario Outline: Intron Retention
		Given an intron retention event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating intron retention is refseq_curated
		When the intron retention events are annotated with the annotation dataset
		Then the resulting annotations of intron retention should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event               | intron | comment |
			| chr9  | 34646787 | 34647088 | +      | ir   | NM_000155.4    | intron 1 retention  | 1      |         |
			| chrX  | 40054043 | 40054255 | -      | ir   | NM_001123385.2 | intron 13 retention | 13     |         |


	Scenario Outline: Cryptic Acceptors
		Given a cryptic acceptor event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating cryptic acceptors is refseq_curated
		When the cryptic acceptor events are annotated with the annotation dataset
		Then the resulting annotations of cryptic acceptors should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event                              | intron | comment  |
			| chr9  | 34647259 | 34647474 | +      | sj   | NM_000155.4    | cryptic intron 2 acceptor @ -18    | 2      | intronic |
			| chr9  | 34647259 | 34647474 | +      | sj   | NM_001258332.2 | cryptic intron 2 acceptor @ -358   | 2      | intronic |
			| chrX  | 41144627 | 41146305 | +      | sj   | NM_001039591.3 | cryptic intron 11 acceptor @ -2064 | 11     | intronic |
			| chr9  | 34646787 | 34647120 | +      | sj   | NM_000155.4    | cryptic exon 2 acceptor @ +32      | 1      | exonic   |
			| chrX  | 41144627 | 41148373 | +      | sj   | NM_001039591.3 | cryptic exon 12 acceptor @ +5      | 11     | exonic   |
			| chrX  | 41223403 | 41224806 | +      | sj   | NM_001039591.3 | cryptic exon 40 acceptor @ +65     | 39     | exonic   |

	# event_type | region type | region number | distance_from_annotated

	Scenario Outline: Cryptic Donors
		Given a cryptic donor event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating cryptic donors is refseq_curated
		When the cryptic donor events are annotated with the annotation dataset
		Then the resulting annotations of cryptic donors should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event                           | intron | comment                                                                                                                    |
			| chr9  | 34646920 | 34647088 | +      | sj   | NM_000155.4    | cryptic intron 1 donor @ +134   | 1      | intronic                                                                                                                   |
			| chrX  | 41114762 | 41123470 | +      | sj   | NM_001039591.3 | cryptic intron 1 donor @ +28653 | 1      | intronic                                                                                                                   |
			| chrX  | 41141436 | 41143290 | +      | sj   | NM_001039591.3 | cryptic intron 9 donor @ +5     | 9      | intronic; super interesting, uses an AT acceptor at -4 but soft clipped to +4 donor position, not a canonical donor site 1 |
			| chrX  | 41144347 | 41144521 | +      | sj   | NM_001039591.3 | cryptic intron 10 donor @ +904  | 10     | intronic                                                                                                                   |
			| chrX  | 41216059 | 41217219 | +      | sj   | NM_001039591.3 | cryptic exon 35 donor @ -594    | 35     | exonic; it's definitely -594, I don't know why it fails                                                                                                                     |


	Scenario Outline: Skip-Cryp Events
		Given a skip-cryp event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating skip-cryps is refseq_curated
		When the skip-cryp events are annotated with the annotation dataset
		Then the resulting annotations of skip-cryps should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event                                               | intron     | comment                                   |
			| chrX  | 41148576 | 41152932 | +      | sj   | NM_001039591.3 | exon 13 skipping/cryptic intron 13 acceptor @ -16   | 12, 13     | single skipping/cryptic intronic acceptor |
			| chr9  | 34647259 | 34648091 | +      | sj   | NM_000155.4    | exon 3-4-5 skipping/cryptic intron 5 acceptor @ -24 | 2, 3, 4, 5 | triple skipping/cryptic intronic acceptor |
			| chr9  | 34647259 | 34648047 | +      | sj   | NM_000155.4    | exon 3-4-5 skipping/cryptic intron 5 acceptor @ -68 | 2, 3, 4, 5 | triple skipping/cryptic intronic acceptor |
			| chrX  | 41198618 | 41205302 | +      | sj   | NM_001039591.3 | exon 31 skipping/cryptic exon 30 donor @ -133       | 30, 31     | single skipping/cryptic exonic donor      |


	Scenario Outline: Unannotated Junctions
		Given an unannotated junction with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating unannotated junctions is refseq_curated
		When the unannotated junctions are annotated with the annotation dataset
		Then the resulting annotations of unannotated junctions should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event                                                    | intron | comment                                                             |
			| chrX  | 41149755 | 41151932 | +      | sj   | NM_001039591.3 | unannotated intronic junction                            |        | intronic; this one also skips exon 13                               |
			| chrX  | 41101571 | 41102912 | -      | sj   | NM_001039591.3 | unannotated intronic junction (opposite strand)          | 1      | intronic; opposite strand                                           |
			| chrX  | 40081196 | 40081258 | +      | sj   | NM_001123385.2 | unannotated intronic junction (opposite strand)          | 1      | intronic; opposite strand                                           |
			| chrX  | 41229636 | 41230503 | +      | sj   | NM_001039591.3 | unannotated exonic junction                              |        | exonic                                                              |
			| chrX  | 40096604 | 40097680 | +      | sj   | NM_001123385.2 | unannotated intronic/exonic junction (opposite strand)   | 1      | intronic/exonic; opposite strand                                    |
			| chrX  | 40097722 | 40098095 | +      | sj   | NM_001123385.2 | unannotated intergenic/exonic junction (opposite strand) | NA     | intergenic/exonic; opposite strand                                  |
			| chrX  | 40113978 | 40208244 | -      | sj   | NM_001123385.2 | unannotated intergenic junction                          | NA     | unknown; within alternate intron 1 - & intergenic (mis-map warning) |
			| chrX  | 40115468 | 40117489 | +      | sj   | NM_001123385.2 | unannotated intergenic junction                          | NA     | unknown; within alternate intron 1 - (genuine splicing?)            |


	Scenario Outline: Annotate events using RefSeq MANE
		Given a splicing event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the transcript annotation is refseq_mane
		When the splicing events are annotated
		Then the result should be event <event> for transcript <transcript>

		Examples:
			| test | chrom | start    | end      | strand | type | transcript     | event                                                  | intron | comment                      | comment2                                                                                                         |
			| 1    | chr9  | 34646787 | 34647088 | +      | sj   | NM_000155.4    | canonical exon 1-2 splicing                            | 1      | galt_intron_1_merged_sjir_4  | canonical +                                                                                                      |
			| 2    | chr9  | 34646787 | 34647088 | +      | ir   | NM_000155.4    | intron 1 retention                                     | 1      | galt_intron_1_merged_sjir_2  | intron retention +                                                                                               |
			| 3    | chr9  | 34647259 | 34647831 | +      | sj   | NM_000155.4    | exon 3-4 skipping                                      | 2      | galt_intron_1_merged_sjir_9  | double exon skipping +                                                                                           |
			| 4    | chr9  | 34647259 | 34648114 | +      | sj   | NM_000155.4    | exon 3-4-5 skipping                                    | 2      | galt_intron_1_merged_sjir_6  | triple exon skipping +                                                                                           |
			| 5    | chr9  | 34646920 | 34647088 | +      | sj   | NM_000155.4    | intronic cryptic donor @ +134                          | 1      | galt_intron_1_merged_sjir_0  | cryptic donor exonic +                                                                                           |
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
			| 18   | chrX  | 40077970 | 40177006 | -      | sj   | NM_001123385.2 | intergenic cryptic donor @ -79792                      | 1      | galt_chrX_subset_BCOR        | alternate donor exon 1-2 splicing (NM_001123383.1) (alternate donor) -                                           |
			| 19   | chrX  | 41086110 | 41123470 | +      | sj   | NM_001039591.3 | canonical exon 1-2 splicing                            | 1      | galt_chrX_subset_USP9X       | canonical +                                                                                                      |
			| 20   | chrX  | 41087666 | 41123470 | +      | sj   | NM_001039591.3 | intronic cryptic donor @ +1557                         | 1      | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 21   | chrX  | 41101571 | 41102912 | -      | sj   | NM_001039591.3 | unannotated intronic junction (opposite strand)        | 1      | galt_chrX_subset_USP9X       | 2                                                                                                                |
			| 22   | chrX  | 41114762 | 41123470 | +      | sj   | NM_001039591.3 | intronic cryptic donor @ +28653                        | 1      | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 23   | chrX  | 41123725 | 41128999 | +      | sj   | NM_001039591.3 | canonical exon 2-3 splicing                            | 2      | galt_chrX_subset_USP9X       | 2                                                                                                                |
			| 24   | chrX  | 41134838 | 41140655 | +      | sj   | NM_001039591.3 | exon 6 skipping                                        |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 25   | chrX  | 41141436 | 41143290 | +      | sj   | NM_001039591.3 | intronic cryptic donor @ +5                            |        | galt_chrX_subset_USP9X       | super interesting, uses an AT acceptor at -4 but soft clipped to +4 donor position, not a canonical donor site 1 |
			| 26   | chrX  | 41144347 | 41144521 | +      | sj   | NM_001039591.3 | intronic cryptic donor @ +904                          |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 27   | chrX  | 41144627 | 41146305 | +      | sj   | NM_001039591.3 | intronic cryptic acceptor @ -2064                      |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 28   | chrX  | 41144627 | 41148373 | +      | sj   | NM_001039591.3 | exonic cryptic acceptor @ +5                           |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 29   | chrX  | 41148576 | 41152932 | +      | sj   | NM_001039591.3 | exon 13 skipping/intronic cryptic acceptor @ -16       |        | galt_chrX_subset_USP9X       | cryptic annotated in other transcripts but not with skipping                                                     |
			| 30   | chrX  | 41149755 | 41151932 | +      | sj   | NM_001039591.3 | unannotated intronic junction                          |        | galt_chrX_subset_USP9X       | this one also skips exon 13                                                                                      |
			| 31   | chrX  | 41153082 | 41167481 | +      | sj   | NM_001039591.3 | exon 15-16 skipping                                    |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 32   | chrX  | 41198618 | 41205302 | +      | sj   | NM_001039591.3 | exon 31 skipping/exonic cryptic donor @ -133           |        | galt_chrX_subset_USP9X       | actually an exonic cryptic donor plus skipping of exon 31                                                        |
			| 33   | chrX  | 41216059 | 41217219 | +      | sj   | NM_001039591.3 | exonic cryptic donor @ -594                            |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 34   | chrX  | 41223403 | 41224806 | +      | sj   | NM_001039591.3 | exonic cryptic acceptor @ +65                          |        | galt_chrX_subset_USP9X       | 1                                                                                                                |
			| 35   | chrX  | 41229636 | 41230503 | +      | sj   | NM_001039591.3 | unannotated exonic junction                            |        | galt_chrX_subset_USP9X       | 0                                                                                                                |