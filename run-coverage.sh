#!/bin/bash

coverage run --source=./src --branch -m behave --tags=~@wip
coverage report
coverage lcov -o coverage.xml
coverage html