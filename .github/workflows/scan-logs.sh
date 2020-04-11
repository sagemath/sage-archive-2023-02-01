#!/usr/bin/env bash
# The markup in the output is a GitHub Actions logging command
# https://help.github.com/en/actions/automating-your-workflow-with-github-actions/development-tools-for-github-actions
LOGS=${1-logs}
find "$LOGS" -type f -name "*.log" -exec sh -c '\
    if tail -100 "{}" 2>/dev/null | grep "^Error" >/dev/null; then \
        echo :":"error file={}:":" ==== ERROR IN LOG FILE {} ====; cat "{}" ; \
    elif tail -20 "{}" 2>/dev/null | grep -E "^(Warning: Error testing|make: .*test.*Error)" >/dev/null; then \
        echo :":"warning file={}:":" ==== TESTSUITE FAILURE IN LOG FILE {} ====; cat "{}"; \
    fi' \;
