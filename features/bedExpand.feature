Feature: Expand collapsed annotation in UCSC Transcript Bed File

    Scenario: RefSeq MANE transcript bed is parsed into its expanded form
        Given the RefSeq MANE transcript bed file is
            | chr  | start  | end    | name                                     | score | strand |
            | chr1 | 65433  | 65519  | NM_001005484.2_intron_0_0_chr1_65434_f   | 0     | +      |
            | chr1 | 65573  | 69036  | NM_001005484.2_intron_1_0_chr1_65574_f   | 0     | +      |
            | chr1 | 924948 | 925921 | NM_001385641.1_intron_0_0_chr1_924949_f  | 0     | +      |
            | chr1 | 926013 | 930154 | NM_001385641.1_intron_1_0_chr1_926014_f  | 0     | +      |
            | chr1 | 930336 | 931038 | NM_001385641.1_intron_2_0_chr1_930337_f  | 0     | +      |
            | chr1 | 931089 | 935771 | NM_001385641.1_intron_3_0_chr1_931090_f  | 0     | +      |
            | chr1 | 935896 | 939039 | NM_001385641.1_intron_4_0_chr1_935897_f  | 0     | +      |
            | chr1 | 939129 | 939274 | NM_001385641.1_intron_5_0_chr1_939130_f  | 0     | +      |
            | chr1 | 939412 | 941143 | NM_001385641.1_intron_6_0_chr1_939413_f  | 0     | +      |
            | chr1 | 941306 | 942135 | NM_001385641.1_intron_7_0_chr1_941307_f  | 0     | +      |
            | chr1 | 942251 | 942409 | NM_001385641.1_intron_8_0_chr1_942252_f  | 0     | +      |
            | chr1 | 942488 | 942558 | NM_001385641.1_intron_9_0_chr1_942489_f  | 0     | +      |
            | chr1 | 943058 | 943252 | NM_001385641.1_intron_10_0_chr1_943059_f | 0     | +      |
            | chr1 | 943377 | 943697 | NM_001385641.1_intron_11_0_chr1_943378_f | 0     | +      |
            | chr1 | 943808 | 943907 | NM_001385641.1_intron_12_0_chr1_943809_f | 0     | +      |
        When the RefSeq MANE transcript bed file is parsed
        Then the RefSeq MANE transcript bed file should be
            | chr | start  | end    | transcript     | intron | strand |
            | chr1  | 65433  | 65519  | NM_001005484.2 | 1      | +      |
            | chr1  | 65573  | 69036  | NM_001005484.2 | 2      | +      |
            | chr1  | 924948 | 925921 | NM_001385641.1 | 1      | +      |
            | chr1  | 926013 | 930154 | NM_001385641.1 | 2      | +      |
            | chr1  | 930336 | 931038 | NM_001385641.1 | 3      | +      |
            | chr1  | 931089 | 935771 | NM_001385641.1 | 4      | +      |
            | chr1  | 935896 | 939039 | NM_001385641.1 | 5      | +      |
            | chr1  | 939129 | 939274 | NM_001385641.1 | 6      | +      |
            | chr1  | 939412 | 941143 | NM_001385641.1 | 7      | +      |
            | chr1  | 941306 | 942135 | NM_001385641.1 | 8      | +      |
            | chr1  | 942251 | 942409 | NM_001385641.1 | 9      | +      |
            | chr1  | 942488 | 942558 | NM_001385641.1 | 10     | +      |
            | chr1  | 943058 | 943252 | NM_001385641.1 | 11     | +      |
            | chr1  | 943377 | 943697 | NM_001385641.1 | 12     | +      |
            | chr1  | 943808 | 943907 | NM_001385641.1 | 13     | +      |