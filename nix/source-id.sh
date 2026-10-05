# shellcheck shell=bash
# The sources the build compiles. With jj, which snapshots the working copy
# first: an empty @'s one parent, a commit that outlives the run, or else
# @'s own commit id, which covers the uncommitted edits but is rewritten by
# the next edit, and dirty: true. With git alone, HEAD, or the commit
# `git stash create` makes of the tracked edits on top of it (touching
# nothing), and dirty: true. Prints the id, then dirty.
source_id() {
  local stash
  if jj --version >/dev/null 2>&1 && jj root >/dev/null 2>&1; then
    jj log --no-graph -r @ -T 'if(empty && parents.len() == 1,
      parents.map(|c| c.commit_id()).join(""), commit_id)
      ++ " " ++ if(empty, "false", "true") ++ "\n"'
  elif git rev-parse HEAD >/dev/null 2>&1; then
    stash=$(git stash create)
    if [[ -n $stash || -n $(git status --porcelain) ]]; then
      echo "${stash:-$(git rev-parse HEAD)} true"
    else
      echo "$(git rev-parse HEAD) false"
    fi
  else
    echo "unknown false"
  fi
}
