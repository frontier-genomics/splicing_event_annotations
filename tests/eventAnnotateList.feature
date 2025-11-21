Feature: Annotate a list of splicing events using RefSeq Curated Annotations

    Scenario: Annotate list of splicing events
        Given the splicing events
            | event | chrom | start    | end      | strand | type | transcript     |
            | 1     | chr9  | 34646787 | 34647088 | +      | sj   | NM_000155.4    |
            | 2     | chrX  | 40054043 | 40054255 | -      | sj   | NM_001123385.2 |
            | 3     | chrX  | 41086110 | 41123470 | +      | sj   | NM_001039591.3 |
            | 4     | chrX  | 41123725 | 41128999 | +      | sj   | NM_001039591.3 |
            | 5     | chrX  | 13749534 | 13751248 | +      | sj   | NM_003611.3    |
            | 6     | chrX  | 13751369 | 13753367 | +      | sj   | NM_003611.3    |
            | 7     | chrX  | 13749534 | 13751248 | +      | sj   | NM_001330210.2 |
            | 8     | chrX  | 13751369 | 13753367 | +      | sj   | NM_001330210.2 |
            | 9     | chrX  | 13749534 | 13753367 | +      | sj   | NM_001330209.2 |
            | 10    | chrX  | 53256062 | 53291894 | -      | sj   | NM_001111125.3 |
            | 11    | chrX  | 53256062 | 53291894 | -      | sj   | NM_001410736.1 |
            | 12    | chrX  | 53279660 | 53281495 | -      | sj   | NM_001243197.2 |
            | 13    | chrX  | 53267065 | 53279563 | -      | sj   | NM_001243197.2 |
            | 14    | chrX  | 53279660 | 53281495 | -      | sj   | NM_015075.2    |
            | 15    | chrX  | 53256062 | 53279563 | -      | sj   | NM_015075.2    |
            | 16    | chrX  | 53235185 | 53236321 | -      | sj   | NM_001111125.3 |
            | 17    | chrX  | 13749534 | 13753367 | +      | sj   | NM_003611.3    |
            | 18    | chrX  | 13749534 | 13753367 | +      | sj   | NM_001330210.2 |
            | 19    | chr9  | 34647259 | 34647831 | +      | sj   | NM_000155.4    |
            | 20    | chrX  | 41229828 | 41230500 | +      | sj   | NM_001039591.3 |
            | 21    | chrX  | 41229828 | 41230500 | +      | sj   | NM_001410749.1 |
            | 22    | chrX  | 13751369 | 13753367 | +      | sj   | NM_001330209.2 |
            | 23    | chrX  | 53256062 | 53279563 | -      | sj   | NM_001111125.3 |
            | 24    | chrX  | 53256062 | 53279563 | -      | sj   | NM_001410736.1 |
            | 25    | chrX  | 41087666 | 41123470 | +      | sj   | NM_001039591.3 |
            | 26    | chrX  | 13749534 | 13751248 | +      | sj   | NM_001330209.2 |
            | 27    | chrX  | 53267065 | 53279563 | -      | sj   | NM_015075.2    |
            | 28    | chrX  | 53256062 | 53291894 | -      | sj   | NM_015075.2    |
            | 29    | chrX  | 40077970 | 40177006 | -      | sj   | NM_001123385.2 |
            | 30    | chrX  | 53256062 | 53279563 | -      | sj   | NM_001243197.2 |
            | 31    | chrX  | 53279660 | 53281495 | -      | sj   | NM_001111125.3 |
            | 32    | chrX  | 53267065 | 53279563 | -      | sj   | NM_001111125.3 |
            | 33    | chrX  | 53279660 | 53281495 | -      | sj   | NM_001410736.1 |
            | 34    | chrX  | 53267065 | 53279563 | -      | sj   | NM_001410736.1 |
            | 35    | chrX  | 53256062 | 53291894 | -      | sj   | NM_001243197.2 |
            | 36    | chrX  | 40057322 | 40062745 | -      | sj   | NM_001123385.2 |
            | 37    | chrX  | 41134838 | 41140655 | +      | sj   | NM_001039591.3 |
            | 38    | chrX  | 41153082 | 41167481 | +      | sj   | NM_001039591.3 |
            | 39    | chrX  | 53238307 | 53243331 | -      | sj   | NM_001111125.3 |
            | 40    | chr9  | 34647259 | 34648114 | +      | sj   | NM_000155.4    |
            | 41    | chrX  | 53238307 | 53246968 | -      | sj   | NM_001111125.3 |
            | 42    | chr9  | 34646787 | 34647088 | +      | ir   | NM_000155.4    |
            | 43    | chrX  | 40054043 | 40054255 | -      | ir   | NM_001123385.2 |
            | 44    | chr9  | 34647259 | 34647474 | +      | sj   | NM_000155.4    |
            | 45    | chr9  | 34647259 | 34647474 | +      | sj   | NM_001258332.2 |
            | 46    | chrX  | 41144627 | 41146305 | +      | sj   | NM_001039591.3 |
            | 47    | chrX  | 53247226 | 53248113 | -      | sj   | NM_001111125.3 |
            | 48    | chr9  | 34646787 | 34647120 | +      | sj   | NM_000155.4    |
            | 49    | chrX  | 41144627 | 41148373 | +      | sj   | NM_001039591.3 |
            | 50    | chrX  | 41223403 | 41224806 | +      | sj   | NM_001039591.3 |
            | 51    | chrX  | 53247110 | 53248113 | -      | sj   | NM_001111125.3 |
            | 52    | chr9  | 34646920 | 34647088 | +      | sj   | NM_000155.4    |
            | 53    | chrX  | 41114762 | 41123470 | +      | sj   | NM_001039591.3 |
            | 54    | chrX  | 41141436 | 41143290 | +      | sj   | NM_001039591.3 |
            | 55    | chrX  | 41144347 | 41144521 | +      | sj   | NM_001039591.3 |
            | 56    | chrX  | 53247136 | 53248044 | -      | sj   | NM_001111125.3 |
            | 57    | chrX  | 41216059 | 41217219 | +      | sj   | NM_001039591.3 |
            | 58    | chrX  | 41217284 | 41218371 | +      | sj   | NM_001039591.3 |
            | 59    | chrX  | 53247136 | 53248139 | -      | sj   | NM_001111125.3 |
            | 60    | chrX  | 41148576 | 41152932 | +      | sj   | NM_001039591.3 |
            | 61    | chr9  | 34647259 | 34648091 | +      | sj   | NM_000155.4    |
            | 62    | chr9  | 34647259 | 34648047 | +      | sj   | NM_000155.4    |
            | 63    | chrX  | 53243517 | 53248720 | -      | sj   | NM_001111125.3 |
            | 64    | chrX  | 53243419 | 53248720 | -      | sj   | NM_001111125.3 |
            | 65    | chr9  | 34647259 | 34647892 | +      | sj   | NM_000155.4    |
            | 66    | chrX  | 53243472 | 53248044 | -      | sj   | NM_001111125.3 |
            | 67    | chr9  | 34647601 | 34647831 | +      | sj   | NM_000155.4    |
            | 68    | chrX  | 53243472 | 53248139 | -      | sj   | NM_001111125.3 |
            | 69    | chrX  | 41198618 | 41205302 | +      | sj   | NM_001039591.3 |
            | 70    | chrX  | 53235185 | 53236474 | -      | sj   | NM_001111125.3 |
            | 71    | chrX  | 41149755 | 41151932 | +      | sj   | NM_001039591.3 |
            | 72    | chrX  | 41101571 | 41102912 | -      | sj   | NM_001039591.3 |
            | 73    | chrX  | 40081196 | 40081258 | +      | sj   | NM_001123385.2 |
            | 74    | chrX  | 41229636 | 41230503 | +      | sj   | NM_001039591.3 |
            | 75    | chrX  | 40096604 | 40097680 | +      | sj   | NM_001123385.2 |
            | 76    | chrX  | 40097722 | 40098095 | +      | sj   | NM_001123385.2 |
            | 77    | chrX  | 40096723 | 40098095 | +      | sj   | NM_001123385.2 |
            | 78    | chrX  | 40113978 | 40208244 | -      | sj   | NM_001123385.2 |
            | 79    | chrX  | 40115468 | 40117489 | +      | sj   | NM_001123385.2 |
        And the annotation dataset is refseq
        When the list of splicing events is annotated
        Then the resulting annotations should be
            | event | event_type                           | location            | distance_from_authentic | event                                                                | intron       | comment                                                                                                          | transcript     |
            | 1     | canonical                            | NA                  | NA                      | canonical exon 1-2 splicing                                          | 1            | canonical splicing                                                                                               | NM_000155.4    |
            | 2     | canonical                            | NA                  | NA                      | canonical exon 13-14 splicing                                        | 13           | canonical splicing                                                                                               | NM_001123385.2 |
            | 3     | canonical                            | NA                  | NA                      | canonical exon 1-2 splicing                                          | 1            | canonical splicing                                                                                               | NM_001039591.3 |
            | 4     | canonical                            | NA                  | NA                      | canonical exon 2-3 splicing                                          | 2            | canonical splicing                                                                                               | NM_001039591.3 |
            | 5     | canonical                            | NA                  | NA                      | canonical exon 9-10 splicing                                         | 9            | canonical splicing                                                                                               | NM_003611.3    |
            | 6     | canonical                            | NA                  | NA                      | canonical exon 10-11 splicing                                        | 10           | canonical splicing                                                                                               | NM_003611.3    |
            | 7     | canonical                            | NA                  | NA                      | canonical exon 10-11 splicing                                        | 10           | canonical splicing                                                                                               | NM_001330210.2 |
            | 8     | canonical                            | NA                  | NA                      | canonical exon 11-12 splicing                                        | 11           | canonical splicing                                                                                               | NM_001330210.2 |
            | 9     | canonical                            | NA                  | NA                      | canonical exon 9-10 splicing                                         | 9            | canonical splicing                                                                                               | NM_001330209.2 |
            | 10    | canonical                            | NA                  | NA                      | canonical exon 2-3 splicing                                          | 2            | canonical splicing                                                                                               | NM_001111125.3 |
            | 11    | canonical                            | NA                  | NA                      | canonical exon 2-3 splicing                                          | 2            | canonical splicing                                                                                               | NM_001410736.1 |
            | 12    | canonical                            | NA                  | NA                      | canonical exon 1-2 splicing                                          | 1            | canonical splicing                                                                                               | NM_001243197.2 |
            | 13    | canonical                            | NA                  | NA                      | canonical exon 2-3 splicing                                          | 2            | canonical splicing                                                                                               | NM_001243197.2 |
            | 14    | canonical                            | NA                  | NA                      | canonical exon 1-2 splicing                                          | 1            | canonical splicing                                                                                               | NM_015075.2    |
            | 15    | canonical                            | NA                  | NA                      | canonical exon 2-3 splicing                                          | 2            | canonical splicing                                                                                               | NM_015075.2    |
            | 16    | alternate exon skipping              | NA                  | NA                      | alternate exon 14 skipping (NM_001410736.1, NM_015075.2)             | 13, 14       | single exon skipping; 2 alt transcripts                                                                          | NM_001111125.3 |
            | 17    | alternate exon skipping              | NA                  | NA                      | alternate exon 10 skipping (NM_001330209.2 exon 9-10)                | 9, 10        | single exon skipping; 1 alt transcript                                                                           | NM_003611.3    |
            | 18    | alternate exon skipping              | NA                  | NA                      | alternate exon 11 skipping (NM_001330209.2 exon 9-10)                | 10, 11       | single exon skipping; 1 alt transcript                                                                           | NM_001330210.2 |
            | 19    | alternate exon skipping              | NA                  | NA                      | alternate exon 3-4 skipping (NM_001258332.2 exon 2-3)                | 2, 3, 4      | double exon skipping; 1 alt transcript                                                                           | NM_000155.4    |
            | 20    | alternate donor                      | intron 43           | 49                      | alternate intron 43 donor @ +49 (NM_001039590.3, NM_001410748.1)     | 43           | intronic donor; 2 alt transcripts                                                                                | NM_001039591.3 |
            | 21    | alternate donor                      | intron 44           | 49                      | alternate intron 44 donor @ +49 (NM_001039590.3, NM_001410748.1)     | 44           | intronic donor; 2 alt transcripts                                                                                | NM_001410749.1 |
            | 22    | alternate donor                      | intron 9            | 1836                    | alternate intron 9 donor @ +1836 (NM_001330210.2, NM_003611.3)       | 9            | intronic donor; 2 alt transcripts                                                                                | NM_001330209.2 |
            | 23    | alternate donor                      | intron 2            | 12332                   | alternate intron 2 donor @ +12332 (NM_015075.2 exon 2-3)             | 2            | intronic donor; 1 alt transcript                                                                                 | NM_001111125.3 |
            | 24    | alternate donor                      | intron 2            | 12332                   | alternate intron 2 donor @ +12332 (NM_015075.2 exon 2-3)             | 2            | intronic donor; 1 alt transcript                                                                                 | NM_001410736.1 |
            | 25    | alternate donor                      | intron 1            | 1557                    | alternate intron 1 donor @ +1557 (NM_001410748.1, NM_001410749.1)    | 1            | intronic donor; 2 alt transcripts                                                                                | NM_001039591.3 |
            | 26    | alternate acceptor                   | intron 9            | -2120                   | alternate intron 9 acceptor @ -2120 (NM_001330210.2, NM_003611.3)    | 9            | intronic acceptor; 2 alt transcripts                                                                             | NM_001330209.2 |
            | 27    | alternate acceptor                   | intron 2            | -11004                  | alternate intron 2 acceptor @ -11004 (NM_001243197.2 exon 2-3)       | 2            | intronic acceptor; 1 alt transcript                                                                              | NM_015075.2    |
            | 28    | alternate donor                      | intergenic          | -12331                  | alternate intergenic donor @ -12331 (NM_001111125.3, NM_001410736.1) | NA           | intergenic donor; 2 alt transcripts                                                                              | NM_015075.2    |
            | 29    | alternate donor                      | intergenic          | -79792                  | alternate intergenic donor @ -79792 (NM_001123383.1, NM_001123384.2) | NA           | intergenic donor; 2 alt transcripts                                                                              | NM_001123385.2 |
            | 30    | alternate acceptor                   | intergenic          | 11003                   | alternate intergenic acceptor @ +11003 (NM_015075.2 exon 2-3)        | NA           | intergenic acceptor; 1 alt transcript                                                                            | NM_001243197.2 |
            | 31    | alternate intron                     | intronic            | NA                      | alternate intronic junction (NM_001243197.2, NM_015075.2)            | NA           | unannotated intronic junction; 2 alt transcripts                                                                 | NM_001111125.3 |
            | 32    | alternate intron                     | intronic            | NA                      | alternate intronic junction (NM_001243197.2 exon 2-3)                | NA           | unannotated intronic junction; 1 alt transcript                                                                  | NM_001111125.3 |
            | 33    | alternate intron                     | intronic            | NA                      | alternate intronic junction (NM_001243197.2, NM_015075.2)            | NA           | unannotated intronic junction; 2 alt transcripts                                                                 | NM_001410736.1 |
            | 34    | alternate intron                     | intronic            | NA                      | alternate intronic junction (NM_001243197.2 exon 2-3)                | NA           | unannotated intronic junction; 1 alt transcript                                                                  | NM_001410736.1 |
            | 35    | alternate intron                     | intergenic          | NA                      | alternate intergenic junction (NM_001111125.3, NM_001410736.1)       | NA           | unannotated intergenic junction; 2 alt transcripts                                                               | NM_001243197.2 |
            | 36    | exon skipping                        | NA                  | NA                      | exon 10 skipping                                                     | 9, 10        | single skipping                                                                                                  | NM_001123385.2 |
            | 37    | exon skipping                        | NA                  | NA                      | exon 6 skipping                                                      | 5, 6         | single skipping                                                                                                  | NM_001039591.3 |
            | 38    | exon skipping                        | NA                  | NA                      | exon 15-16 skipping                                                  | 14, 15, 16   | double skipping                                                                                                  | NM_001039591.3 |
            | 39    | exon skipping                        | NA                  | NA                      | exon 10-11 skipping                                                  | 9, 10, 11    | double skipping                                                                                                  | NM_001111125.3 |
            | 40    | exon skipping                        | NA                  | NA                      | exon 3-4-5 skipping                                                  | 2, 3, 4, 5   | triple skipping                                                                                                  | NM_000155.4    |
            | 41    | exon skipping                        | NA                  | NA                      | exon 9-10-11 skipping                                                | 8, 9, 10, 11 | triple skipping                                                                                                  | NM_001111125.3 |
            | 42    | intron retention                     | NA                  | NA                      | intron 1 retention                                                   | 1            | intron retention                                                                                                 | NM_000155.4    |
            | 43    | intron retention                     | NA                  | NA                      | intron 13 retention                                                  | 13           | intron retention                                                                                                 | NM_001123385.2 |
            | 44    | cryptic acceptor                     | intron 2            | -18                     | cryptic intron 2 acceptor @ -18                                      | 2            | cryptic intronic acceptor                                                                                        | NM_000155.4    |
            | 45    | cryptic acceptor                     | intron 2            | -358                    | cryptic intron 2 acceptor @ -358                                     | 2            | cryptic intronic acceptor                                                                                        | NM_001258332.2 |
            | 46    | cryptic acceptor                     | intron 11           | -2064                   | cryptic intron 11 acceptor @ -2064                                   | 11           | cryptic intronic acceptor                                                                                        | NM_001039591.3 |
            | 47    | cryptic acceptor                     | intron 7            | -91                     | cryptic intron 7 acceptor @ -91                                      | 7            | cryptic intronic acceptor                                                                                        | NM_001111125.3 |
            | 48    | cryptic acceptor                     | exon 2              | 32                      | cryptic exon 2 acceptor @ +32                                        | 1            | cryptic exonic acceptor                                                                                          | NM_000155.4    |
            | 49    | cryptic acceptor                     | exon 12             | 5                       | cryptic exon 12 acceptor @ +5                                        | 11           | cryptic exonic acceptor                                                                                          | NM_001039591.3 |
            | 50    | cryptic acceptor                     | exon 40             | 65                      | cryptic exon 40 acceptor @ +65                                       | 39           | cryptic exonic acceptor                                                                                          | NM_001039591.3 |
            | 51    | cryptic acceptor                     | exon 8              | 26                      | cryptic exon 8 acceptor @ +26                                        | 7            | cryptic exonic acceptor                                                                                          | NM_001111125.3 |
            | 52    | cryptic donor                        | intron 1            | 134                     | cryptic intron 1 donor @ +134                                        | 1            | cryptic intronic donor                                                                                           | NM_000155.4    |
            | 53    | cryptic donor                        | intron 1            | 28653                   | cryptic intron 1 donor @ +28653                                      | 1            | cryptic intronic donor                                                                                           | NM_001039591.3 |
            | 54    | cryptic donor                        | intron 9            | 5                       | cryptic intron 9 donor @ +5                                          | 9            | super interesting, uses an AT acceptor at -4 but soft clipped to +4 donor position, not a canonical donor site 1 | NM_001039591.3 |
            | 55    | cryptic donor                        | intron 10           | 904                     | cryptic intron 10 donor @ +904                                       | 10           | cryptic intronic donor                                                                                           | NM_001039591.3 |
            | 56    | cryptic donor                        | intron 7            | 70                      | cryptic intron 7 donor @ +70                                         | 7            | cryptic intronic donor                                                                                           | NM_001111125.3 |
            | 57    | cryptic donor                        | exon 35             | -594                    | cryptic exon 35 donor @ -594                                         | 35           | cryptic exonic donor                                                                                             | NM_001039591.3 |
            | 58    | cryptic donor                        | exon 36             | -60                     | cryptic exon 36 donor @ -60                                          | 36           | cryptic exonic donor                                                                                             | NM_001039591.3 |
            | 59    | cryptic donor                        | exon 7              | -26                     | cryptic exon 7 donor @ -26                                           | 7            | cryptic exonic donor                                                                                             | NM_001111125.3 |
            | 60    | exon skipping, cryptic acceptor      | intron 13           | -16                     | exon 13 skipping/cryptic intron 13 acceptor @ -16                    | 12, 13       | single skipping/cryptic intronic acceptor                                                                        | NM_001039591.3 |
            | 61    | exon skipping, cryptic acceptor      | intron 5            | -24                     | exon 3-4-5 skipping/cryptic intron 5 acceptor @ -24                  | 2, 3, 4, 5   | triple skipping/cryptic intronic acceptor                                                                        | NM_000155.4    |
            | 62    | exon skipping, cryptic acceptor      | intron 5            | -68                     | exon 3-4-5 skipping/cryptic intron 5 acceptor @ -68                  | 2, 3, 4, 5   | triple skipping/cryptic intronic acceptor                                                                        | NM_000155.4    |
            | 63    | exon skipping, cryptic acceptor      | intron 8            | -46                     | exon 7-8 skipping/cryptic intron 8 acceptor @ -46                    | 6, 7, 8      | double skipping/cryptic intronic acceptor                                                                        | NM_001111125.3 |
            | 64    | exon skipping, cryptic acceptor      | exon 9              | 53                      | exon 7-8 skipping/cryptic exon 9 acceptor @ +53                      | 6, 7, 8      | double skipping/cryptic exonic acceptor                                                                          | NM_001111125.3 |
            | 65    | exon skipping, cryptic acceptor      | exon 5              | 61                      | exon 3-4 skipping/cryptic exon 5 acceptor @ +61                      | 2, 3, 4      | double skipping/cryptic exonic acceptor                                                                          | NM_000155.4    |
            | 66    | exon skipping, cryptic donor         | intron 7            | 70                      | exon 8 skipping/cryptic intron 7 donor @ +70                         | 7, 8         | single skipping/cryptic intronic donor                                                                           | NM_001111125.3 |
            | 67    | exon skipping, cryptic donor         | intron 3            | 34                      | exon 4 skipping/cryptic intron 3 donor @ +34                         | 3, 4         | single skipping/cryptic intronic donor                                                                           | NM_000155.4    |
            | 68    | exon skipping, cryptic donor         | exon 7              | -26                     | exon 8 skipping/cryptic exon 7 donor @ -26                           | 7, 8         | single skipping/cryptic exonic donor                                                                             | NM_001111125.3 |
            | 69    | exon skipping, cryptic donor         | exon 30             | -133                    | exon 31 skipping/cryptic exon 30 donor @ -133                        | 30, 31       | single skipping/cryptic exonic donor                                                                             | NM_001039591.3 |
            | 70    | exon skipping, cryptic donor         | exon 13             | -153                    | exon 14 skipping/cryptic exon 13 donor @ -153                        | 13, 14       | single skipping/cryptic exonic donor                                                                             | NM_001111125.3 |
            | 71    | unannotated intron                   | intronic            | NA                      | unannotated intronic junction                                        | NA           | intronic; this one also skips exon 13                                                                            | NM_001039591.3 |
            | 72    | unannotated intron (opposite strand) | intronic            | NA                      | unannotated intronic junction (opposite strand)                      | NA           | intronic; opposite strand                                                                                        | NM_001039591.3 |
            | 73    | unannotated intron (opposite strand) | intronic            | NA                      | unannotated intronic junction (opposite strand)                      | NA           | intronic; opposite strand                                                                                        | NM_001123385.2 |
            | 74    | unannotated intron                   | exonic              | NA                      | unannotated exonic junction                                          | NA           | exonic                                                                                                           | NM_001039591.3 |
            | 75    | unannotated intron (opposite strand) | intronic/exonic     | NA                      | unannotated intronic/exonic junction (opposite strand)               | NA           | intronic/exonic; opposite strand                                                                                 | NM_001123385.2 |
            | 76    | unannotated intron (opposite strand) | intergenic/exonic   | NA                      | unannotated intergenic/exonic junction (opposite strand)             | NA           | intergenic/exonic; opposite strand                                                                               | NM_001123385.2 |
            | 77    | unannotated intron (opposite strand) | intergenic/intronic | NA                      | unannotated intergenic/intronic junction (opposite strand)           | NA           | intergenic/intronic; opposite strand                                                                             | NM_001123385.2 |
            | 78    | unannotated intron                   | intergenic          | NA                      | unannotated intergenic junction                                      | NA           | unknown; within alternate intron 1 - & intergenic (mis-map warning)                                              | NM_001123385.2 |
            | 79    | unannotated intron                   | intergenic          | NA                      | unannotated intergenic junction                                      | NA           | unknown; within alternate intron 1 - (genuine splicing?)                                                         | NM_001123385.2 |





















