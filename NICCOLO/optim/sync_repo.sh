#!/usr/bin/env bash
# Pull the team's changes, then push this work to the repo. Safe to re-run: it commits only when something changed.
#
#   bash sync_repo.sh "what changed"
#
# PULL: fetch origin, then merge origin/main into my branch (in my own worktree).
# PUSH: copy the whitelisted files into NICCOLO/optim/ of that worktree, commit, push branch niccolo/optim-variants (never force, never main).
# The shared checkout /home/hackathon/hackathon-team-3 hosts another session's uncommitted work: only `git fetch` ever runs there.
# Large files (the RF3 archive and the unzipped folds, ~370 MB) stay on the VM; nothing over 5 MB and no credential is allowed through.
set -euo pipefail
SRC=/home/hackathon/nick-work/optim
REPO=/home/hackathon/hackathon-team-3
WT=/home/hackathon/wt/optim
BR=niccolo/optim-variants
DEST=NICCOLO/optim
MSG=${1:-"optim variants: sync"}

git -C "$REPO" fetch origin --prune --quiet                                    # PULL, part 1
[ -d "$WT" ] || git -C "$REPO" worktree add -q -b "$BR" "$WT" origin/main
cd "$WT"
[ "$(git rev-parse --abbrev-ref HEAD)" = "$BR" ] || { echo "worktree $WT is not on $BR"; exit 1; }
git merge --no-edit --quiet origin/main || { git merge --abort; echo "merge with origin/main conflicts: resolve by hand in $WT"; exit 1; }   # PULL, part 2
echo "base: origin/main $(git rev-parse --short origin/main) merged into $BR"

mkdir -p "$DEST"
rsync -a --ignore-missing-args --files-from=- "$SRC/" "$WT/$DEST/" <<'LIST'
README.md
sync_repo.sh
gen_score.py
build_submission.py
score_structures.py
rescore_current_prefilter.py
rank_final.py
jobs/submit_japanfold.py
jobs/rf3_payload.json
jobs/boltz2_payload.json
jobs/targets_manifest.csv
jobs/rf3_payload.job_id
jobs/results/rf3_payload/job.json
jobs/results/rf3_payload/results.json
out/final_variants.csv
out/final_variants_scored.csv
out/final_variants_rescored_current_prefilter.csv
out/final_variants_ranked.csv
out/final_variants_ranked_abguard.csv
out/rf3_structure_scores.csv
out/selected_diverse_200.csv
out/designs_candidates.csv
out/check_fast.csv
out/check_final.csv
out/check_ranked_fast.csv
out/pool_stats.json
out/pool_top5000_by_pred_kd_before_hole_check.csv
out/gen_score.log
out/score_structures.log
out/rescore_current_prefilter.log
out/rank_final.log
LIST

# guards: nothing large, no credential
big=$(find "$DEST" -type f -size +5M | head -3); [ -z "$big" ] || { echo "refusing: files over 5 MB: $big"; exit 1; }
if [ -n "${JAPANFOLD_API_KEY:-}" ] && grep -rqF -f <(printf '%s\n' "$JAPANFOLD_API_KEY") "$DEST"; then echo "refusing: the JapanFold key is in the files"; exit 1; fi
if grep -rqE 'Bearer [A-Za-z0-9._-]{20,}' "$DEST"; then echo "refusing: a bearer-token pattern is in the files"; exit 1; fi

git add -A "$DEST"
if git diff --cached --quiet; then echo "nothing new to commit"; else
  git commit -q -m "$MSG" -m "Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
  echo "committed: $(git log --oneline -1)"
fi
git push -u origin "$BR"                                                       # PUSH
echo "pushed $BR at $(git rev-parse --short HEAD)"
