Feature: Annotate a tsv of splicing events using RefSeq Curated Annotations
@wip
    Scenario: Annotate tsv of splicing events
        Given the splicing event tsv resources/test_input.tsv
        And the annotation dataset for the tsv is refseq
        When the tsv of splicing events is annotated
        Then the resulting annotated tsv should be resources/test_output.tsv
            












