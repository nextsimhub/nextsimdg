#!/bin/bash
# Bash script for preparing a file to be processed by Doxygen.
#
# In particular, it adds the @file and @date entries to the frontmatter of the file.

set -eu

# Check for valid input
if [ -z "$1" ]; then
  echo "Error: The doxy_prep.sh script requires a filename as input."
  exit 1
fi
FILE="$1"

# Check if the file exists
if [ ! -f "${FILE}" ]; then
  echo "Error: File '${FILE}' does not exist."
  exit 1
fi

# Check if the file is tracked by git
if ! git ls-files --error-unmatch "${FILE}" >/dev/null 2>&1; then
  echo "Error: File '${FILE}' is not tracked by git."
  exit 1
fi

# Only consider files containing '@author'
if grep -q '@author' "${FILE}"; then

  # Get the base filename
  FILE_NO_PATH=$(basename "${FILE}")

  # Get the date of the most recent change using git blame
  LAST_CHANGE_DATE=$(git blame --line-porcelain "${FILE}" | grep "^author-time " | sort -nr | head -n 1 | awk '{print $2}' | xargs -I{} date -d @{} +"%d %b %Y")
  if [ -z "${LAST_CHANGE_DATE}" ]; then
    echo "Error: Unable to determine the last change date for '${FILE}'."
    exit 1
  fi

  # Insert frontmatter @file and @date entries required by Doxygen, separated by a blank line
  sed -i "1a \ * @file    ${FILE_NO_PATH}\n *\n * @date    ${LAST_CHANGE_DATE}" "${FILE}"
fi
