#!/bin/bash

# Run behave but ignore tests that are marked as WIP
behave --tags=~@wip