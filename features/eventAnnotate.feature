Feature: Annotate events in cortar data using RefSeq Curated Annotations

	Scenario Outline: Canonical Splicing
		Given a canonical splicing event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating canonical splicing is refseq
		When the canonical splicing events are annotated with the annotation dataset
		Then the resulting annotations of canonical splicing should be event <event> of type <event_type> at intron <intron>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event_type | event                         | intron | comment |
			| chr9  | 34646787 | 34647088 | +      | sj   | NM_000155.4    | canonical  | canonical exon 1-2 splicing   | 1      |         |
			| chrX  | 40054043 | 40054255 | -      | sj   | NM_001123385.2 | canonical  | canonical exon 13-14 splicing | 13     |         |
			| chrX  | 41086110 | 41123470 | +      | sj   | NM_001039591.3 | canonical  | canonical exon 1-2 splicing   | 1      |         |
			| chrX  | 41123725 | 41128999 | +      | sj   | NM_001039591.3 | canonical  | canonical exon 2-3 splicing   | 2      |         |
			| chrX  | 13749534 | 13751248 | +      | sj   | NM_003611.3    | canonical  | canonical exon 9-10 splicing  | 9      |         |
			| chrX  | 13751369 | 13753367 | +      | sj   | NM_003611.3    | canonical  | canonical exon 10-11 splicing | 10     |         |
			| chrX  | 13749534 | 13751248 | +      | sj   | NM_001330210.2 | canonical  | canonical exon 10-11 splicing | 10     |         |
			| chrX  | 13751369 | 13753367 | +      | sj   | NM_001330210.2 | canonical  | canonical exon 11-12 splicing | 11     |         |
			| chrX  | 13749534 | 13753367 | +      | sj   | NM_001330209.2 | canonical  | canonical exon 9-10 splicing  | 9      |         |
			| chrX  | 53256062 | 53291894 | -      | sj   | NM_001111125.3 | canonical  | canonical exon 2-3 splicing   | 2      |         |
			| chrX  | 53256062 | 53291894 | -      | sj   | NM_001410736.1 | canonical  | canonical exon 2-3 splicing   | 2      |         |
			| chrX  | 53279660 | 53281495 | -      | sj   | NM_001243197.2 | canonical  | canonical exon 1-2 splicing   | 1      |         |
			| chrX  | 53267065 | 53279563 | -      | sj   | NM_001243197.2 | canonical  | canonical exon 2-3 splicing   | 2      |         |
			| chrX  | 53279660 | 53281495 | -      | sj   | NM_015075.2    | canonical  | canonical exon 1-2 splicing   | 1      |         |
			| chrX  | 53256062 | 53279563 | -      | sj   | NM_015075.2    | canonical  | canonical exon 2-3 splicing   | 2      |         |


	Scenario Outline: Alternate Splicing
		Given an alternate splicing event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating alternate splicing is refseq
		When the alternate splicing events are annotated with the annotation dataset
		Then the resulting annotations of alternate splicing should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event_type              | location   | distance_from_authentic | event                                                                | intron | comment                                            |
			| chrX  | 53235185 | 53236321 | -      | sj   | NM_001111125.3 | alternate exon skipping |            |                         | alternate exon 14 skipping (NM_001410736.1, NM_015075.2)             |        | single exon skipping; 2 alt transcripts            |
			| chrX  | 13749534 | 13753367 | +      | sj   | NM_003611.3    | alternate exon skipping |            |                         | alternate exon 10 skipping (NM_001330209.2 exon 9-10)                |        | single exon skipping; 1 alt transcript             |
			| chrX  | 13749534 | 13753367 | +      | sj   | NM_001330210.2 | alternate exon skipping |            |                         | alternate exon 11 skipping (NM_001330209.2 exon 9-10)                |        | single exon skipping; 1 alt transcript             |
			| chr9  | 34647259 | 34647831 | +      | sj   | NM_000155.4    | alternate exon skipping |            |                         | alternate exon 3-4 skipping (NM_001258332.2 exon 2-3)                | 2      | double exon skipping; 1 alt transcript             |
			| chrX  | 41229828 | 41230500 | +      | sj   | NM_001039591.3 | alternate donor         | intronic   | +49                     | alternate intron 43 donor @ +49 (NM_001039590.3, NM_001410748.1)     | 43     | intronic donor; 2 alt transcripts                  |
			| chrX  | 41229828 | 41230500 | +      | sj   | NM_001410749.1 | alternate donor         | intronic   | +49                     | alternate intron 44 donor @ +49 (NM_001039590.3, NM_001410748.1)     | 44     | intronic donor; 2 alt transcripts                  |
			| chrX  | 13751369 | 13753367 | +      | sj   | NM_001330209.2 | alternate donor         | intronic   | +1836                   | alternate intron 9 donor @ +1836 (NM_001330210.2, NM_003611.3)       | 9      | intronic donor; 2 alt transcripts                  |
			| chrX  | 53256062 | 53279563 | -      | sj   | NM_001111125.3 | alternate donor         | intronic   | +12332                  | alternate intron 2 donor @ +12332 (NM_015075.2 exon 2-3)             | 2      | intronic donor; 1 alt transcript                   |
			| chrX  | 53256062 | 53279563 | -      | sj   | NM_001410736.1 | alternate donor         | intronic   | +12332                  | alternate intron 2 donor @ +12332 (NM_015075.2 exon 2-3)             |        | intronic donor; 1 alt transcript                   |
			| chrX  | 41087666 | 41123470 | +      | sj   | NM_001039591.3 | alternate donor         | intronic   | +1557                   | alternate intron 1 donor @ +1557 (NM_001410748.1, NM_001410749.1)    | 1      | intronic donor; 2 alt transcripts                  |
			| chrX  | 13749534 | 13751248 | +      | sj   | NM_001330209.2 | alternate acceptor      | intronic   | -2120                   | alternate intron 9 acceptor @ -2120 (NM_001330210.2, NM_003611.3)    | 9      | intronic acceptor; 2 alt transcripts               |
			| chrX  | 53267065 | 53279563 | -      | sj   | NM_015075.2    | alternate acceptor      | intronic   | -11004                  | alternate intron 2 acceptor @ -11004 (NM_001243197.2 exon 2-3)       |        | intronic acceptor; 1 alt transcript                |
			| chrX  | 53256062 | 53291894 | -      | sj   | NM_015075.2    | alternate donor         | intergenic | -12331                  | alternate intergenic donor @ -12331 (NM_001111125.3, NM_001410736.1) |        | intergenic donor; 2 alt transcripts                |
			| chrX  | 40077970 | 40177006 | -      | sj   | NM_001123385.2 | alternate donor         | intergenic | -79792                  | alternate intergenic donor @ -79792 (NM_001123383.1, NM_001123384.2) | 1      | intergenic donor; 2 alt transcripts                |
			| chrX  | 53256062 | 53279563 | -      | sj   | NM_001243197.2 | alternate acceptor      | intergenic | +11003                  | alternate intergenic acceptor @ +11003 (NM_015075.2 exon 2-3)        |        | intergenic acceptor; 1 alt transcript              |
			| chrX  | 53279660 | 53281495 | -      | sj   | NM_001111125.3 | alternate intron        | intronic   |                         | alternate intronic junction (NM_001243197.2, NM_015075.2)            |        | unannotated intronic junction; 2 alt transcripts   |
			| chrX  | 53267065 | 53279563 | -      | sj   | NM_001111125.3 | alternate intron        | intronic   |                         | alternate intronic junction (NM_001243197.2 exon 2-3)                |        | unannotated intronic junction; 1 alt transcript    |
			| chrX  | 53279660 | 53281495 | -      | sj   | NM_001410736.1 | alternate intron        | intronic   |                         | alternate intronic junction (NM_001243197.2, NM_015075.2)            |        | unannotated intronic junction; 2 alt transcripts   |
			| chrX  | 53267065 | 53279563 | -      | sj   | NM_001410736.1 | alternate intron        | intronic   |                         | alternate intronic junction (NM_001243197.2 exon 2-3)                |        | unannotated intronic junction; 1 alt transcript    |
			| chrX  | 53256062 | 53291894 | -      | sj   | NM_001243197.2 | alternate intron        | intergenic |                         | alternate intergenic junction (NM_001111125.3, NM_001410736.1)       |        | unannotated intergenic junction; 2 alt transcripts |


	Scenario Outline: Exon Skipping
		Given an exon skipping event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating exon skipping is refseq
		When the exon skipping events are annotated with the annotation dataset
		Then the resulting annotations of exon skipping should be event <event> of type <event_type> at intron <intron>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event_type    | event                 | intron       | comment |
			| chrX  | 40057322 | 40062745 | -      | sj   | NM_001123385.2 | exon skipping | exon 10 skipping      | 9, 10        | single  |
			| chrX  | 41134838 | 41140655 | +      | sj   | NM_001039591.3 | exon skipping | exon 6 skipping       | 5, 6         | single  |
			| chrX  | 41153082 | 41167481 | +      | sj   | NM_001039591.3 | exon skipping | exon 15-16 skipping   | 14, 15, 16   | double  |
			| chrX  | 53238307 | 53243331 | -      | sj   | NM_001111125.3 | exon skipping | exon 10-11 skipping   | 9, 10, 11    | double  |
			| chr9  | 34647259 | 34648114 | +      | sj   | NM_000155.4    | exon skipping | exon 3-4-5 skipping   | 2, 3, 4, 5   | triple  |
			| chrX  | 53238307 | 53246968 | -      | sj   | NM_001111125.3 | exon skipping | exon 9-10-11 skipping | 8, 9, 10, 11 | triple  |


	Scenario Outline: Intron Retention
		Given an intron retention event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating intron retention is refseq
		When the intron retention events are annotated with the annotation dataset
		Then the resulting annotations of intron retention should be event <event> of type <event_type> at intron <intron>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event_type       | event               | intron | comment |
			| chr9  | 34646787 | 34647088 | +      | ir   | NM_000155.4    | intron retention | intron 1 retention  | 1      |         |
			| chrX  | 40054043 | 40054255 | -      | ir   | NM_001123385.2 | intron retention | intron 13 retention | 13     |         |


	Scenario Outline: Cryptic Acceptors
		Given a cryptic acceptor event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating cryptic acceptors is refseq
		When the cryptic acceptor events are annotated with the annotation dataset
		Then the resulting annotations of cryptic acceptors should be event <event> of location <location> and type <event_type> at intron <intron> at a distance of <distance_from_authentic>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event_type       | location | distance_from_authentic | event                              | intron | comment |
			| chr9  | 34647259 | 34647474 | +      | sj   | NM_000155.4    | cryptic acceptor | intronic | -18                     | cryptic intron 2 acceptor @ -18    | 2      |         |
			| chr9  | 34647259 | 34647474 | +      | sj   | NM_001258332.2 | cryptic acceptor | intronic | -358                    | cryptic intron 2 acceptor @ -358   | 2      |         |
			| chrX  | 41144627 | 41146305 | +      | sj   | NM_001039591.3 | cryptic acceptor | intronic | -2064                   | cryptic intron 11 acceptor @ -2064 | 11     |         |
			| chrX  | 53247226 | 53248113 | -      | sj   | NM_001111125.3 | cryptic acceptor | intronic | -91                     | cryptic intron 7 acceptor @ -91    | 7      |         |
			| chr9  | 34646787 | 34647120 | +      | sj   | NM_000155.4    | cryptic acceptor | exonic   | +32                     | cryptic exon 2 acceptor @ +32      | 1      |         |
			| chrX  | 41144627 | 41148373 | +      | sj   | NM_001039591.3 | cryptic acceptor | exonic   | +5                      | cryptic exon 12 acceptor @ +5      | 11     |         |
			| chrX  | 41223403 | 41224806 | +      | sj   | NM_001039591.3 | cryptic acceptor | exonic   | +65                     | cryptic exon 40 acceptor @ +65     | 39     |         |
			| chrX  | 53247110 | 53248113 | -      | sj   | NM_001111125.3 | cryptic acceptor | exonic   | +26                     | cryptic exon 8 acceptor @ +26      | 7      |         |


	Scenario Outline: Cryptic Donors
		Given a cryptic donor event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating cryptic donors is refseq
		When the cryptic donor events are annotated with the annotation dataset
		Then the resulting annotations of cryptic donors should be event <event> of location <location> and type <event_type> at intron <intron> at a distance of <distance_from_authentic>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event_type    | location | distance_from_authentic | event                           | intron | comment                                                                                                          |
			| chr9  | 34646920 | 34647088 | +      | sj   | NM_000155.4    | cryptic donor | intronic | +134                    | cryptic intron 1 donor @ +134   | 1      |                                                                                                                  |
			| chrX  | 41114762 | 41123470 | +      | sj   | NM_001039591.3 | cryptic donor | intronic | +28653                  | cryptic intron 1 donor @ +28653 | 1      |                                                                                                                  |
			| chrX  | 41141436 | 41143290 | +      | sj   | NM_001039591.3 | cryptic donor | intronic | +5                      | cryptic intron 9 donor @ +5     | 9      | super interesting, uses an AT acceptor at -4 but soft clipped to +4 donor position, not a canonical donor site 1 |
			| chrX  | 41144347 | 41144521 | +      | sj   | NM_001039591.3 | cryptic donor | intronic | +904                    | cryptic intron 10 donor @ +904  | 10     |                                                                                                                  |
			| chrX  | 53247136 | 53248044 | -      | sj   | NM_001111125.3 | cryptic donor | intronic | +70                     | cryptic intron 7 donor @ +70    | 7      |                                                                                                                  |
			| chrX  | 41216059 | 41217219 | +      | sj   | NM_001039591.3 | cryptic donor | exonic   | -594                    | cryptic exon 35 donor @ -594    | 35     |                                                                                                                  |
			| chrX  | 41217284 | 41218371 | +      | sj   | NM_001039591.3 | cryptic donor | exonic   | -60                     | cryptic exon 36 donor @ -60     | 36     |                                                                                                                  |
			| chrX  | 53247136 | 53248139 | -      | sj   | NM_001111125.3 | cryptic donor | exonic   | -26                     | cryptic exon 7 donor @ -26      | 7      |                                                                                                                  |

	Scenario Outline: Skip-Cryp Events
		Given a skip-cryp event with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating skip-cryps is refseq
		When the skip-cryp events are annotated with the annotation dataset
		Then the resulting annotations of skip-cryps should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event_type                      | location | distance_from_authentic | event                                               | intron     | comment                                   |
			| chrX  | 41148576 | 41152932 | +      | sj   | NM_001039591.3 | exon skipping, cryptic acceptor | intronic | -16                     | exon 13 skipping/cryptic intron 13 acceptor @ -16   | 12, 13     | single skipping/cryptic intronic acceptor |
			| chr9  | 34647259 | 34648091 | +      | sj   | NM_000155.4    | exon skipping, cryptic acceptor | intronic | -24                     | exon 3-4-5 skipping/cryptic intron 5 acceptor @ -24 | 2, 3, 4, 5 | triple skipping/cryptic intronic acceptor |
			| chr9  | 34647259 | 34648047 | +      | sj   | NM_000155.4    | exon skipping, cryptic acceptor | intronic | -68                     | exon 3-4-5 skipping/cryptic intron 5 acceptor @ -68 | 2, 3, 4, 5 | triple skipping/cryptic intronic acceptor |
			| chrX  | 53243517 | 53248720 | -      | sj   | NM_001111125.3 | exon skipping, cryptic acceptor | intronic | -46                     | exon 7-8 skipping/cryptic intron 8 acceptor @ -46   | 6, 7, 8    | double skipping/cryptic intronic acceptor |
			| chrX  | 53243419 | 53248720 | -      | sj   | NM_001111125.3 | exon skipping, cryptic acceptor | exonic   | +53                     | exon 7-8 skipping/cryptic exon 9 acceptor @ +53     | 6, 7, 8    | double skipping/cryptic exonic acceptor   |
			| chr9  | 34647259 | 34647892 | +      | sj   | NM_000155.4    | exon skipping, cryptic acceptor | exonic   | +61                     | exon 3-4 skipping/cryptic exon 5 acceptor @ +61     | 2, 3, 4    | double skipping/cryptic exonic acceptor   |
			| chrX  | 53243472 | 53248044 | -      | sj   | NM_001111125.3 | exon skipping, cryptic donor    | intronic | +70                     | exon 8 skipping/cryptic intron 7 donor @ +70        | 7, 8, 9    | single skipping/cryptic intronic donor    |
			| chr9  | 34647601 | 34647831 | +      | sj   | NM_000155.4    | exon skipping, cryptic donor    | intronic | +34                     | exon 4 skipping/cryptic intron 3 donor @ +34        | 7, 8, 9    | single skipping/cryptic intronic donor    |
			| chrX  | 53243472 | 53248139 | -      | sj   | NM_001111125.3 | exon skipping, cryptic donor    | exonic   | -26                     | exon 8 skipping/cryptic exon 7 donor @ -26          | 7, 8, 9    | single skipping/cryptic exonic donor      |
			| chrX  | 41198618 | 41205302 | +      | sj   | NM_001039591.3 | exon skipping, cryptic donor    | exonic   | -133                    | exon 31 skipping/cryptic exon 30 donor @ -133       | 30, 31     | single skipping/cryptic exonic donor      |
			| chrX  | 53235185 | 53236474 | -      | sj   | NM_001111125.3 | exon skipping, cryptic donor    | exonic   | -153                    | exon 14 skipping/cryptic exon 13 donor @ -153       | 13, 14     | single skipping/cryptic exonic donor      |


	Scenario Outline: Unannotated Junctions
		Given an unannotated junction with chrom <chrom>, start <start>, end <end>, strand <strand>, transcript <transcript>, and type <type>
		And the annotation dataset to be used for annotating unannotated junctions is refseq
		When the unannotated junctions are annotated with the annotation dataset
		Then the resulting annotations of unannotated junctions should be event <event>

		Examples:
			| chrom | start    | end      | strand | type | transcript     | event_type                           | location          | event                                                    | intron | comment                                                             |
			| chrX  | 41149755 | 41151932 | +      | sj   | NM_001039591.3 | unannotated intron                   | intronic          | unannotated intronic junction                            |        | intronic; this one also skips exon 13                               |
			| chrX  | 41101571 | 41102912 | -      | sj   | NM_001039591.3 | unannotated intron (opposite strand) | intronic          | unannotated intronic junction (opposite strand)          | 1      | intronic; opposite strand                                           |
			| chrX  | 40081196 | 40081258 | +      | sj   | NM_001123385.2 | unannotated intron (opposite strand) | intronic          | unannotated intronic junction (opposite strand)          | 1      | intronic; opposite strand                                           |
			| chrX  | 41229636 | 41230503 | +      | sj   | NM_001039591.3 | unannotated intron                   | exonic            | unannotated exonic junction                              |        | exonic                                                              |
			| chrX  | 40096604 | 40097680 | +      | sj   | NM_001123385.2 | unannotated intron (opposite strand) | intronic/exonic   | unannotated intronic/exonic junction (opposite strand)   | 1      | intronic/exonic; opposite strand                                    |
			| chrX  | 40097722 | 40098095 | +      | sj   | NM_001123385.2 | unannotated intron (opposite strand) | intergenic/exonic | unannotated intergenic/exonic junction (opposite strand) | NA     | intergenic/exonic; opposite strand                                  |
			| chrX  | 40113978 | 40208244 | -      | sj   | NM_001123385.2 | unannotated intron                   | intergenic        | unannotated intergenic junction                          | NA     | unknown; within alternate intron 1 - & intergenic (mis-map warning) |
			| chrX  | 40115468 | 40117489 | +      | sj   | NM_001123385.2 | unannotated intron                   | intergenic        | unannotated intergenic junction                          | NA     | unknown; within alternate intron 1 - (genuine splicing?)            |


	Scenario Outline: Fetch overlapping RefSeq MANE transcripts for event
		Given a splicing event with chrom <chrom>, start <start>, end <end>, and strand <strand>
		When the overlapping transcripts are fetched
		Then the result should be transcript <transcript>, gene <gene>, and warning <warning>

		Examples:
			| test | chrom | start    | end      | strand | transcript     | gene    | warning                                               |
			| 1    | chr9  | 34646787 | 34647088 | +      | NM_000155.4    | GALT    | none                                                  |
			| 2    | chr9  | 34646787 | 34647088 | +      | NM_000155.4    | GALT    | none                                                  |
			| 3    | chr9  | 34647259 | 34647831 | +      | NM_000155.4    | GALT    | none                                                  |
			| 4    | chr9  | 34647259 | 34648114 | +      | NM_000155.4    | GALT    | none                                                  |
			| 5    | chr9  | 34646920 | 34647088 | +      | NM_000155.4    | GALT    | none                                                  |
			| 6    | chr9  | 34646787 | 34647120 | +      | NM_000155.4    | GALT    | none                                                  |
			| 7    | chr9  | 34647259 | 34647474 | +      | NM_000155.4    | GALT    | none                                                  |
			| 8    | chr9  | 34647259 | 34648091 | +      | NM_000155.4    | GALT    | none                                                  |
			| 9    | chr9  | 34647259 | 34648047 | +      | NM_000155.4    | GALT    | none                                                  |
			| 10   | chrX  | 40054043 | 40054255 | -      | NM_001123385.2 | BCOR    | none                                                  |
			| 11   | chrX  | 40054043 | 40054255 | -      | NM_001123385.2 | BCOR    | none                                                  |
			| 12   | chrX  | 40057322 | 40062745 | -      | NM_001123385.2 | BCOR    | none                                                  |
			| 13   | chrX  | 40081196 | 40081258 | +      | NM_001123385.2 | BCOR    | MANE transcript found on opposite strand only.        |
			| 14   | chrX  | 40096604 | 40097680 | +      | NM_001123385.2 | BCOR    | MANE transcript found on opposite strand only.        |
			| 15   | chrX  | 40097722 | 40098095 | +      | NM_001123385.2 | BCOR    | MANE transcript found on opposite strand only.        |
			| 16   | chrX  | 40113978 | 40208244 | -      | unknown        | unknown | no MANE transcript match found                        |
			| 17   | chrX  | 40115468 | 40117489 | +      | unknown        | unknown | no MANE transcript match found                        |
			| 18   | chrX  | 40077970 | 40177006 | -      | NM_001123385.2 | BCOR    | none                                                  |
			| 19   | chrX  | 41086110 | 41123470 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 20   | chrX  | 41087666 | 41123470 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 21   | chrX  | 41101571 | 41102912 | -      | NM_001039591.3 | USP9X   | MANE transcript found on opposite strand only.        |
			| 22   | chrX  | 41114762 | 41123470 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 23   | chrX  | 41123725 | 41128999 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 24   | chrX  | 41134838 | 41140655 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 25   | chrX  | 41141436 | 41143290 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 26   | chrX  | 41144347 | 41144521 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 27   | chrX  | 41144627 | 41146305 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 28   | chrX  | 41144627 | 41148373 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 29   | chrX  | 41148576 | 41152932 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 30   | chrX  | 41149755 | 41151932 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 31   | chrX  | 41153082 | 41167481 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 32   | chrX  | 41198618 | 41205302 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 33   | chrX  | 41216059 | 41217219 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 34   | chrX  | 41223403 | 41224806 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 35   | chrX  | 41229636 | 41230503 | +      | NM_001039591.3 | USP9X   | none                                                  |
			| 36   | chrX  | 53267065 | 53279563 | +      | NM_001111125.3 | IQSEC2  | MANE transcript found on opposite strand only.        |
			| 37   | chrX  | 54445540 | 54446203 | -      | NM_004463.3    | FGD1    | Overlapping MANE transcript on opposite strand found. |
			| 38   | chrX  | 54448698 | 54449551 | +      | NM_004463.3    | FGD1    | MANE transcript found on opposite strand only.        |