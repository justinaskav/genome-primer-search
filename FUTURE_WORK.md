# Future Work

Items identified during the pipeline audit that are not yet implemented.
Each entry is scoped enough to act on without reading the conversation history
that produced it.

## Current state (so future-you doesn't re-introduce the things we fixed)

Already landed:

- **Error strategy is narrow.** `retry → ignore` only applies to `DOWNLOAD_GENOME`. Every other process fails fast. Failed downloads are written to `reports/failed_genomes.tsv` (only when any fail). Do not re-introduce a global `errorStrategy`.
- **No silent defaults for `--primers`.** It's a required parameter; `main.nf` errors out with an actionable message if missing or pointing at a nonexistent file. Same for `--genomes`.
- **Filter parser carries structure.** `bin/filter_amplicons.py` emits a `stats.json` with `totals`, `per_primer`, `per_contig`, and a full `amplicons[]` array. Per-amplicon records include forward/reverse primer hit objects with `primer_seq`, `strand`, `pos`, `mismatches`, plus `sequence_id` (contig), `sequence_desc`, and `filter_status`. The HTML report consumes this directly — no more proportional allocation.
- **`max_amplicon_size` default is 10 kb.** Aligned with the README. The filter explanation no longer uses the "circular artifact" framing; it describes cross-locus pairings, which is what's actually happening.
- **Thermo Tm/ΔTm/P(amp) are decoupled from references.** `bin/thermodynamic_analysis.py` always computes these from the primer sequence and master-mix conditions. The reference FASTA's only job is target-vs-off-target classification, and only when supplied.
- **`is_intended_target` is tri-state.** `True` = reference declared and binding region matches (≥95 % identity via local alignment); `False` = reference declared and rejected; `None` = no reference for this primer. The HTML report has an explicit intent banner and hides specificity / off-target sections when no references are provided. The specificity TSV writes a header-only file in that case.

---

## T0.4 — Replace the specificity score

**Files:** `bin/generate_thermo_report.py`

The current per-primer "Specificity Score" is `100 * (1 - weighted_offtargets / total_hits)`
with weights `high=1.0, medium=0.5, low=0.1` (around `generate_thermo_report.py:200`).
The weights are arbitrary, the bins (Excellent/Good/Fair/Poor) are arbitrary, and
the score conflates off-target count with severity — a primer with 100 low-risk
off-targets scores similarly to one with 10.

**Replace with** a per-primer panel showing the constituent numbers:

- Intended-target hits: count, mean ΔTm, mean P(amp)
- Off-target hits: count, max P(amp), max ΔTm
- 3'-end-mismatch hits (any non-zero in `forward_3prime_mm_est` / `reverse_3prime_mm_est`): count

A "good" primer reads off this table directly: many intended-target hits with
high P(amp), few off-target hits, especially few with high P(amp). No score,
no rating bins.

Only render this panel when intent was declared (`metadata.references_provided`).
Helpers `is_target` / `is_off_target` / `is_unclassified` at the top of
`generate_thermo_report.py` already exist for the filters.

**Verify:** run with `--enable_thermo_analysis --thermo_references <real>.fasta`
against a genome with multiple binding sites for each primer. Confirm:

- The panel renders.
- Numbers add up: `target_hits + offtarget_hits + unclassified_hits == total_hits` per primer.
- No score, no Excellent/Good/Fair/Poor labels remain anywhere in the report or TSVs.

---

## T0.5 — De-cliff the P(amplification) binning

**Files:** `bin/generate_thermo_report.py`

`risk_classification` is currently `high / medium / low` based on hard cutoffs
at 0.7 and 0.3. The HTML renders these as discrete colored chips, so a hit at
P(amp) = 0.699 reads "Low" and a hit at 0.701 reads "High" — same biology,
different label.

**Replace** the chip with a continuous gradient bar inside the cell, e.g.:

```html
<td>
  <span class="p-amp-num">0.687</span>
  <span class="p-amp-bar" style="width: 68.7%; background: linear-gradient(to right, #ffc107, #dc3545);"></span>
</td>
```

Keep the existing color tokens (`#28a745` low, `#ffc107` medium, `#dc3545` high);
just plot the value continuously. The underlying `risk_classification` string
stays in the JSON for downstream consumers.

**Verify:** run against a primer set that produces a spread of P(amp) values
from 0.4 to 0.9. Open the HTML and visually confirm the column is a smooth
gradient, not a discrete chip.

---

## T0.6 — Mark the heuristic dimer/hairpin checks as heuristic

**Files:** `bin/thermodynamic_analysis.py:818–864`, `bin/generate_thermo_report.py`

The dimer scan only checks the 3'-end of the primer against its own
reverse-complement; the hairpin scan flags 3-bp stems regardless of loop
length, so a 50-nt loop with a 3-bp stem reads "high" hairpin risk despite
being thermodynamically unstable.

A full swap to `primer3-py` (`calcHairpin`, `calcHomodimer`, `calcHeterodimer`) is
out of scope for this pass — too much new surface. Minimum-viable cleanup:

1. Tighten `check_hairpin()` to require loop length ≥ 4 nt (the inner `j` loop already starts at `i + 4` but the stem-length scoring downstream doesn't gate on loop reasonableness — tighten the upper bound and reject stems with very long unmatched intervening regions).
2. In `generate_thermo_report.py`, rename the HTML column headings from "Dimer Risk" / "Hairpin Risk" to "Dimer flag (heuristic)" / "Hairpin flag (heuristic)" and add a footnote: *"Heuristic scan based on 3'-end complementarity and palindrome detection. For production primer screening, validate with primer3 or IDT OligoAnalyzer."*

**Verify:** synthetic primer with a known 3-bp stem and 30-nt loop should no
longer be flagged as high hairpin risk. Primer with a known stable hairpin
(e.g. 6-bp stem, 5-nt loop) should still flag.

---

## T0.7 — README rewrite for thermo

**Files:** `README.md` (the "Thermodynamic Off-Target Analysis" section)

The current text claims a perfect-match hit is only classified as an intended
target if a matching reference is provided. That's no longer literally true —
the new logic also requires the binding-region sequence to match the reference
at ≥95 % identity, and the reference's role is no longer to drive Tm/ΔTm. Sync
the docs after T0.4 / T0.5 / T0.6 land.

Touch points:

- "Intended Target vs Off-Target" subsection: reference matches via sequence similarity, not zero-mismatch gating. ΔTm is always defined.
- Behavior with no reference: the report still renders, just without specificity / off-target classification.
- Add a new "What goes in a reference FASTA?" subsection. Entry header must be `>primer_name_reference` (the `_reference` suffix is stripped). Content should be the expected amplicon region as it would appear in the intended-target organism's genome. Tag this with a worked example using `data/references/16S_ref.fasta`.
- After T0.4 lands, drop any mention of the specificity score from the README.

---

## T1.1 — Empty / malformed FASTA sanity check

**Files:** `modules/local/extract_genome.nf:24-26, 39-43`

Current check is `[ ! -s ${input_safe}_genome.fasta ]` — catches zero-byte
files but not files containing only a header, or non-FASTA junk that happens to
be non-empty.

**Fix:** after the existing size check, add a `>` count and fail with a clear
message naming the genome:

```bash
if [ "$(grep -c '^>' ${input_safe}_genome.fasta)" -lt 1 ]; then
    echo "Error: ${input_safe}_genome.fasta contains no FASTA records (${genome_input})" >&2
    exit 1
fi
```

Pairs cleanly with the narrowed `errorStrategy` already in place — `EXTRACT_GENOME`
failures now surface as real errors instead of being retried-then-dropped.

**Verify:** stage an empty FASTA (`echo ">empty" > /tmp/bad.fa`); confirm the
pipeline fails fast at `EXTRACT_GENOME` with the message above and never
reaches `RUN_PRIMERSEARCH`. (You'll need to wire it through a test input that
routes to that file — easiest is to inject a fake download path in a
short-lived test branch.)

---

## T1.2 — Pin conda environments

**Files:** `modules/local/{download_genome,run_primersearch,thermodynamic_analysis}.nf`, new `envs/*.yml`

Current directives read `conda 'bioconda::emboss'` etc. — unpinned. EMBOSS
primersearch output format has been stable, but BioPython's `MeltingTemp`
defaults *have* shifted between minor releases, which silently changes
thermo numbers across runs.

**Fix:**

1. Create `envs/emboss.yml`, `envs/biopython.yml`, `envs/ncbi.yml` with explicit pins:
   ```yaml
   # envs/emboss.yml
   channels: [bioconda, conda-forge]
   dependencies:
     - bioconda::emboss=6.6.0
   ```
2. In each `modules/local/*.nf`, replace the inline conda directive with `conda 'envs/<name>.yml'` (path resolves relative to `projectDir`).
3. Add `bin/requirements.txt`:
   ```
   biopython==1.83
   ```
   (Currently the only pure-Python dependency outside the stdlib.)

Suggested initial pins (cross-check the latest stable in each channel before merging):

| Tool | Pin |
|---|---|
| EMBOSS | `6.6.0` |
| BioPython | `1.83` |
| `ncbi-datasets-cli` | `16.*` |
| `entrez-direct` | `22.*` |

**Verify:** delete the local conda cache (`conda clean --all` or
`rm -rf ~/.conda/envs/<env-name>`), rerun `tests/run_test.sh`. The resolved
versions in `_work/.command.sh` and the conda env manifest must match the
pinned `envs/*.yml`. Multiple runs should produce identical thermo JSON for
the same inputs.

---

## T1.3 — Offline regression test

**Files:** `tests/data/regression_ecoli_16s.fasta` (new), `tests/data/regression_primers.txt` (new), `tests/data/expected_summary.tsv` (new), `tests/test_offline_regression.sh` (new), `tests/run_test.sh` (update)

The repo currently has no end-to-end correctness check — tests download real
genomes from NCBI and verify the pipeline *completes*, but never verify the
numbers. A parser regression or upstream tool drift could change amplicon
counts or sizes without anyone noticing.

**Fix:** bundle a ~50–100 kb fragment of E. coli K-12 (NC_000913.3) around
one 16S rRNA locus (e.g. `rrnB`) as a static FASTA. Run the pipeline against
this fragment with the existing `tests/data/test_primers.txt`. Compare
`_output/reports/summary.tsv` against a checked-in expected file.

Expected: 1 amplicon per primer (V3V4 ≈ 465 bp, 27F ≈ 1506 bp), length within
±1 bp of the published values, zero off-targets in the small fragment.

The test must work offline — feed the FASTA in as a local file (via a small
shim that detects local-path inputs in `DOWNLOAD_GENOME`, or as a separate
test entrypoint that skips that step). Decide which: a local-path branch in
`download_genome.nf` is the more useful primitive but bigger scope; a test-only
skip is faster to land.

Wire into `tests/run_test.sh` so it runs as part of the default sweep.

**Verify:** the test must fail if (a) the EMBOSS version changes its position
arithmetic, (b) the parser drops the mismatch field, (c) the size filter
default changes, or (d) the summary TSV column order changes. The expected
file should be a small, human-readable diff target — not auto-regenerated each
run.

---

## Order of execution

1. **T0.4 → T0.5 → T0.6** in one pass through `bin/generate_thermo_report.py` (and a small dip into `bin/thermodynamic_analysis.py` for the hairpin tightening). Bundle the report changes since they touch the same file.
2. **T0.7** after the above land, so the README describes the final shape rather than a moving target.
3. **T1.3** before **T1.2** — the regression test is the safety net that proves T1.2's version pins didn't change numbers silently. Without T1.3, T1.2 is fail-open.
4. **T1.1** any time — independent of the rest.

T1.2 + T1.3 together are the reproducibility story; together they're worth
their own PR.
