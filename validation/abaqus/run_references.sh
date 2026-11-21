#!/bin/bash
# Abaqus reference results for the validation suite: runs the decks of validation/cases one after the other with
# Abaqus and extracts validation/reference/<case>.json (and <case>.time, wall time of the job in seconds).
# Cases whose name starts with par_ are also run on 4 cpus (<case>__cpus4.json / .time); decks containing the
# comment '** validation: single' also in single precision (<case>__single.json); a deck containing
# '** validation: oldjob=<case>' is run with oldjob=<case> (transfer of results, *IMPORT), after that case; a static
# deck containing '** validation: explicit_static' (elements absent from Abaqus/Standard) is run as an Abaqus/Explicit
# quasi-static analysis converted by static_to_explicit.py, whose final state becomes the reference.
# Run on a machine with an Abaqus licence, then bring validation/reference/ back into the repository.
#
# usage: bash validation/abaqus/run_references.sh [case ...]        (all cases by default)
#        ABAQUS=abq2018 bash validation/abaqus/run_references.sh   (Abaqus command, default: abaqus)
#        FORCE=1 ...                                               (rerun cases that already have a reference)
#        DATACHECK=1 ...                     (input check only: Abaqus pre-processor on every deck, no analysis)
# Explicit decks run in double precision (double=both); work files go to validation/abaqus/work/.
set -u
ABQ=${ABAQUS:-abaqus}
HERE=$(cd "$(dirname "$0")/.." && pwd)
WORK="$HERE/abaqus/work"
mkdir -p "$WORK" "$HERE/reference"
cd "$WORK" || exit 1
cp "$HERE"/cases/*.inc . 2>/dev/null  # included files are read relative to the working directory
cases="$*"
[ -z "$cases" ] && cases=$(ls "$HERE/cases" | sed -n 's/\.inp$//p')
total=$(echo $cases | wc -w)
ok=0
failed=""
i=0

run_one() {  # case, job name, cpus [precision]
  local n=$1 job=$2 cpus=$3 extra="" input="$HERE/cases/$n.inp" qs=""
  : > "$job.run.log"
  if grep -q '^\*\* validation: explicit_static' "$input"; then
    qs=1
    $ABQ python "$HERE/abaqus/static_to_explicit.py" "$input" "$PWD/$job.qs.inp" >> "$job.run.log" 2>&1 || return 1
    input="$PWD/$job.qs.inp"
  fi
  grep -qi "dynamic, explicit" "$input" && extra="double=${4:-both}"
  local old=$(sed -n 's/^\*\* validation: oldjob=//p' "$HERE/cases/$n.inp")
  [ -n "$old" ] && extra="$extra oldjob=$old"
  rm -f "$job.lck"
  local t0=$(date +%s)
  $ABQ job="$job" input="$input" cpus=$cpus $extra interactive ask_delete=OFF >> "$job.run.log" 2>&1
  local rc=$?
  local t1=$(date +%s)
  if [ $rc -eq 0 ] && $ABQ python "$HERE/abaqus/extract_odb.py" "$job.odb" "$HERE/reference/$job.json" \
       "$HERE/cases/manifest.json" "$n" >> "$job.run.log" 2>&1 && { [ -z "$qs" ] || \
       $ABQ python "$HERE/abaqus/static_to_explicit.py" --reduce "$HERE/reference/$job.json" "$HERE/cases/$n.inp" \
       >> "$job.run.log" 2>&1; }; then
    echo "$((t1 - t0))" > "$HERE/reference/$job.time"
    return 0
  fi
  return 1
}

if [ -n "${DATACHECK:-}" ]; then  # syntax check of every deck with the Abaqus pre-processor
  for n in $cases; do
    i=$((i + 1)); rm -f "dc_$n.lck"
    if grep -q '^\*\* validation: oldjob=' "$HERE/cases/$n.inp"; then ok=$((ok + 1)); continue; fi
    input="$HERE/cases/$n.inp"
    if grep -q '^\*\* validation: explicit_static' "$input"; then
      $ABQ python "$HERE/abaqus/static_to_explicit.py" "$input" "$PWD/dc_$n.qs.inp" > /dev/null 2>&1
      input="$PWD/dc_$n.qs.inp"
    fi
    if $ABQ job="dc_$n" input="$input" datacheck interactive ask_delete=OFF > "dc_$n.run.log" 2>&1; then
      ok=$((ok + 1))
    else
      echo "[$i/$total] ERROR  $n: $(cat "dc_$n.dat" "dc_$n.sta" 2>/dev/null | grep -m1 -A2 '[*][*][*]ERROR' | tr -s ' ' | tr '\n' ' ')"
      failed="$failed $n"
    fi
  done
  echo "datacheck: $ok/$total ok"
  [ -z "$failed" ]
  exit
fi
echo "$(date '+%F %T')  $total cases, Abaqus command: $ABQ"
for n in $cases; do
  i=$((i + 1))
  if [ -z "${FORCE:-}" ] && [ -f "$HERE/reference/$n.json" ]; then
    echo "$(date '+%T')  [$i/$total] skip    $n (reference exists)"
    ok=$((ok + 1))
    continue
  fi
  if run_one "$n" "$n" 1; then
    echo "$(date '+%T')  [$i/$total] ok      $n ($(cat "$HERE/reference/$n.time") s)"
    ok=$((ok + 1))
  else
    echo "$(date '+%T')  [$i/$total] FAILED  $n (see validation/abaqus/work/$n.run.log, $n.dat, $n.msg)"
    failed="$failed $n"
  fi
  if grep -q '^\*\* validation: single' "$HERE/cases/$n.inp"; then
    if run_one "$n" "${n}__single" 1 off; then
      echo "$(date '+%T')           ok      ${n}__single ($(cat "$HERE/reference/${n}__single.time") s)"
    else
      echo "$(date '+%T')           FAILED  ${n}__single"
    fi
  fi
  case "$n" in
    par_*)
      if run_one "$n" "${n}__cpus4" 4; then
        echo "$(date '+%T')           ok      ${n}__cpus4 ($(cat "$HERE/reference/${n}__cpus4.time") s)"
      else
        echo "$(date '+%T')           FAILED  ${n}__cpus4"
      fi ;;
  esac
done
echo "$(date '+%F %T')  done: $ok/$total ok"
[ -n "$failed" ] && echo "failed:$failed"
[ -z "$failed" ]
