#!/bin/bash

pytest tests/ --cov=splicing_event_annotator --cov-branch --cov-report=term --cov-report=html --cov-report=lcov:coverage.xml -m "not wip"