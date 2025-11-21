import pytest
from pytest_bdd import given, parsers


@pytest.fixture
def datatable(request):
    """
    Fixture to capture data tables from Gherkin scenarios.
    This will be populated by pytest-bdd when a step has a data table.
    """
    return []


@given(parsers.parse(""), target_fixture="datatable")
def capture_datatable(step):
    """
    Capture the datatable from a step if it exists.
    pytest-bdd passes the step object which may contain a datatable.
    """
    if hasattr(step, 'table'):
        # Convert table to list of dicts for easier access
        table_data = []
        for row in step.table[1:]:  # Skip header
            row_dict = {}
            for i, header in enumerate(step.table[0]):
                row_dict[header] = row[i] if i < len(row) else ''
            table_data.append(row_dict)
        return table_data
    return []
