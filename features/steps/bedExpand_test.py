from behave import given, when, then
from src.bedExpandProcessor import BedExpand
import pandas as pd

@given(u'the RefSeq MANE transcript bed file is')
def step_impl(context):
        context.bed_file = BedExpand(context.table.rows, context.table.headings)

@when(u'the RefSeq MANE transcript bed file is parsed')
def step_impl(context):
    context.expanded = context.bed_file.parse_bed_file()

@then(u'the RefSeq MANE transcript bed file should be')
def step_impl(context):
    context.expanded = pd.DataFrame(context.expanded).astype(str)
    print(context.expanded)
    context.test_data = pd.DataFrame(context.table).astype(str)
    print(context.test_data)
    assert context.expanded.equals(context.test_data)