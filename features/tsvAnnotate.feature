Feature: Annotate a tsv with events
@wip
    Scenario: Annotate test.tsv with events
        Given a tsv resources/test.tsv
        And the annotation dataset is refseq
        When the tsv is annotated
        Then the resulting annotations should be
            | seqnames | start     | end       | strand | gene   | transcript     | intron | event                                                   |
            | chr19    | 38586003  | 38586090  | +      | RYR1   | NM_000540.3    | 103    | canonical exon 103-104 splicing                         |
            | chr4     | 89835693  | 89837097  | -      | SNCA   | NM_001375285.1 | 2      | cryptic donor ~ exon 3                                  |
            | chr19    | 38575962  | 38577917  | +      | RYR1   | NM_000540.3    | 97     | intron 97 retention                                     |
            | chr7     | 150952854 | 150955395 | -      | KCNH2  | NM_000238.4    | 5      | cryptic donor ~ exon 6                                  |
            | chr16    | 2050487   | 2053280   | +      | TSC2   | NM_000548.5    | 3      | exon 3 ~ cryptic acceptor                               |
            | chr9     | 132925740 | 132928766 | -      | TSC1   | NM_000368.5    | 4      | exon 4 skipping                                         |
            | chr11    | 108321421 | 108327644 | +      | ATM    | NM_001351834.2 | 46     | exon 47-48 skipping                                     |
            | chr9     | 132935096 | 132944546 | -      | TSC1   | NM_000368.5    | 1      | cryptic donor ~ exon 2                                  |
            | chr11    | 108225575 | 108227594 | +      | ATM    | NM_001351834.2 | 2      | cryptic donor ~ exon 3                                  |
            | chr9     | 132924907 | 132925586 | -      | TSC1   | NM_000368.5    | 5      | exon 5 ~ cryptic acceptor                               |
            | chr11    | 108293478 | 108299713 | +      | ATM    | NM_001351834.2 | 32     | exon 33-34 skipping                                     |
            | chr12    | 51907696  | 51915224  | +      | ACVRL1 | NM_000020.3    | 3      | exon 2-3-4-5-6 skipping                                 |
            | chr14    | 73146127  | 73147794  | +      | PSEN1  | NM_000021.4    | 1      | cryptic donor ~ exon 2                                  |
            | chr11    | 108229324 | 108243952 | +      | ATM    | NM_001351834.2 | 5      | exon 6 skipping                                         |
            | chr15    | 48441847  | 48472550  | -      | FBN1   | NM_000138.5    | 41     | exon 36-37-38-39-40-41-42-43-44-45-46-47-48-49 skipping |
            | chr21    | 25705259  | 26021846  | -      | APP    | NM_000484.4    | 8      | unannotated junctions                                   |
            | chr21    | 25705259  | 26133639  | -      | APP    | NM_000484.4    | 8      | unannotated junctions                                   |
            | chr2     | 237325096 | 237350141 | -      | COL6A3 | NM_004369.4    | 42     | unannotated junctions                                   |
            | chr21    | 45987659  | 45992176  | +      | COL6A1 | NM_001848.3    | 15     | unannotated junctions                                   |
            | chr2     | 237315245 | 237378665 | -      | COL6A3 | NM_004369.4    | 26     | unannotated junctions                                   |
            | chr2     | 237336213 | 237357811 | -      | COL6A3 | NM_004369.4    | 26     | unannotated junctions                                   |
            | chr21    | 45994215  | 46125265  | +      | COL6A2 | NM_001849.4    | 11     | cryptic donor ~ exon 24                                 |
            | chr21    | 45987176  | 45987609  | +      | COL6A1 | NM_001848.3    | 7      | exon 7 skipping                                         |
            | chr21    | 45975521  | 46112397  | +      | COL6A1 | NM_001848.3    | 27     | unannotated junctions                                   |
            | chr21    | 45997755  | 45998397  | +      | COL6A1 | NM_001848.3    | 22     | cryptic donor ~ exon 24                                 |
            | chr2     | 237315245 | 237378665 | -      | COL6A3 | NM_004369.4    | 14     | unannotated junctions                                   |
            | chr21    | 45984418  | 45986631  | +      | COL6A1 | NM_001848.3    | 3      | unannotated junctions                                   |
            | chr21    | 45992399  | 45998397  | +      | COL6A1 | NM_001848.3    | 20     | exon 19-20-21-22-23 skipping                            |
            | chr15    | 48446823  | 48481654  | -      | FBN1   | NM_000138.5    | 36     | exon 33-34-35-36-37-38-39-40-41-42-43-44-45-46 skipping |




























