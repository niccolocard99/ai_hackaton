#!/usr/bin/env bash
# Cached AMASS searches for the trastuzumab optimisation literature review.
#
# Usage:  AMASS_API_KEY=amass_... bash amass_search.sh [limit]
#
# Every query is cached under cache/<core>__<slug>.json and is NEVER re-fetched
# if the file already exists, so re-running this script costs nothing. Delete a
# cache file to refresh that one query.
#
# Cost: $0.05 per 20 results per search. Default limit 20 => $0.05/query.
# Rate limit: 60 requests / 60 s. We sleep 1.1 s between calls.
set -uo pipefail

: "${AMASS_API_KEY:?set AMASS_API_KEY in the environment (never commit it)}"
BASE=https://api.amass.tech/api/v1
LIMIT=${1:-20}
CACHE="$(dirname "$0")/cache"
mkdir -p "$CACHE"

# core|query
QUERIES=(
  # --- binding affinity / affinity maturation ---
  "biomedcore|trastuzumab affinity maturation CDR mutations"
  "biomedcore|anti-HER2 antibody affinity improvement SPR KD variant"
  "biomedcore|trastuzumab CDRH3 engineering HER2 binding"
  "biomedcore|antibody CDR3 randomization library affinity HER2"
  # --- epitope preservation ---
  "biomedcore|trastuzumab HER2 extracellular domain IV epitope structure"
  "biomedcore|trastuzumab Fab HER2 crystal structure contact residues 1N8Z"
  "biomedcore|HER2 epitope mapping competitive binding trastuzumab pertuzumab"
  # --- thermostability ---
  "biomedcore|anti-HER2 antibody thermostability melting temperature Tm"
  "biomedcore|scFv thermostability engineering framework mutations"
  "biomedcore|antibody affinity stability tradeoff CDR mutations"
  "biomedcore|trastuzumab Fab thermal stability DSC DSF nanoDSF"
  "biomedcore|single chain variable fragment stabilizing mutations disulfide"
  # --- prediction tools / developability ---
  "biomedcore|computational prediction antibody thermostability validation Fv"
  "biomedcore|Aggrescan3D SoluProt NetSolP protein solubility prediction"
  "biomedcore|ProteinMPNN inverse folding antibody design soluble"
  "biomedcore|cell-free protein synthesis antibody expression yield determinants"
  # --- CDRH3 loop structure rules ---
  "biomedcore|CDRH3 loop length conformation canonical structure antibody"
  "biomedcore|de novo antibody loop design RFdiffusion CDR"
  # --- prior art ---
  "patentcore|trastuzumab variant CDR HER2 antibody sequence"
  "patentcore|anti-HER2 antibody affinity matured humanized claim"
  # --- clinical / regulatory context for user-fit framing ---
  "trialcore|trastuzumab HER2 breast cancer resistance"

  # --- BATCH 2: specific structural / biochemical / thermodynamic gaps ---
  # Added after batch 1 showed abstracts carried no trastuzumab thermodynamics,
  # no fragment Tm baseline, and no alanine-scan data.
  "biomedcore|trastuzumab HER2 isothermal titration calorimetry enthalpy entropy binding"
  "biomedcore|trastuzumab alanine scanning mutagenesis HER2 binding hot spot residues"
  "biomedcore|trastuzumab Fab melting temperature nanoDSF domain unfolding biosimilar"
  "biomedcore|trastuzumab scFv single chain construct thermal stability expression"
  "biomedcore|trastuzumab association dissociation rate constants kon koff surface plasmon resonance"
  "biomedcore|HER2 antibody interface buried surface area shape complementarity crystal structure"
  "biomedcore|CDR grafting epitope retention sequence identity threshold antibody"
  "biomedcore|antibody binding thermodynamics enthalpy entropy compensation affinity maturation"
  "patentcore|machine learning computational designed antibody CDR sequence"
)

slug() { echo "$1" | tr '[:upper:] ' '[:lower:]-' | tr -cd 'a-z0-9-' | cut -c1-60; }

fetched=0; cached=0; failed=0
for entry in "${QUERIES[@]}"; do
  core="${entry%%|*}"; query="${entry#*|}"
  out="$CACHE/${core}__$(slug "$query").json"
  if [[ -s "$out" ]]; then
    cached=$((cached+1)); continue
  fi
  code=$(curl -sS --get "$BASE/cores/$core/records" \
    --data-urlencode "query=$query" \
    --data-urlencode "limit=$LIMIT" \
    -H "Authorization: Bearer $AMASS_API_KEY" \
    -w '%{http_code}' -o "$out.tmp" 2>"$out.err")
  if [[ "$code" == "200" ]]; then
    mv "$out.tmp" "$out"; rm -f "$out.err"; fetched=$((fetched+1))
    printf 'ok   %-11s %s\n' "$core" "$query"
  else
    failed=$((failed+1))
    printf 'FAIL %-11s %s -> HTTP %s: %s\n' "$core" "$query" "$code" "$(head -c 200 "$out.tmp" 2>/dev/null)"
    rm -f "$out.tmp"
  fi
  sleep 1.1
done

echo
echo "fetched=$fetched cached=$cached failed=$failed  (cost this run: ~\$$(awk "BEGIN{printf \"%.2f\", $fetched*0.05*$LIMIT/20}"))"
