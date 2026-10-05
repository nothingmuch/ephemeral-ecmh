#!/usr/bin/env bash
# nix run .#validate-commits -- <revision-range>, e.g. main..HEAD
#
# The flake checks at each commit of the range, oldest first. Merge commits
# are skipped: in this history each has the tree of its second parent.

set -euo pipefail

usage() {
  echo "Usage: validate-commits <revision-range>"
}

if [[ ${1:-} == --help ]]; then
  usage
  exit 0
fi

if [[ $# -ne 1 ]]; then
  usage >&2
  exit 2
fi

repo_root=$(git rev-parse --show-toplevel)
commits=$(git rev-list --reverse --no-merges "$1")

if [[ -z $commits ]]; then
  echo "No commits in range."
  exit 0
fi

while read -r commit; do
  echo "Validating $commit"
  nix flake check --no-update-lock-file --print-build-logs \
    "git+file://$repo_root?rev=$commit"
done <<<"$commits"
