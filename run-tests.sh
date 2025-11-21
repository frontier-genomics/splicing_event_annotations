#!/bin/bash

# Run pytest but ignore tests that are marked as WIP
pytest tests/ -m "not wip"