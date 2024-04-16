#!/bin/bash

coverage run --source=./src/python --branch -m behave --tags=~@wip
coverage report
coverage lcov -o coverage.xml
coverage html