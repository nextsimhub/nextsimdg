#!/bin/bash
echo "Generating CHANGELOG" &&
github_changelog_generator -u nextsimdg -p nextsimdg autochangelog-hook -t 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 &&
echo "Running python script" &&
git tag -l -n9 > a.txt &&
python3 clog.py &&
rm a.txt &&
echo "Done! Check CHANGELOG.md"
