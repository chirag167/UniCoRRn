# UnicoRRn

A universal **R**eaction **R**ate converter. UnicoRRn defines X, a generic reaction-rate XML format built on [xmlcoll](https://xmlcoll.readthedocs.io), and converts a lab's instance of it (e.g. X(TAMU)) into other formats used by nuclear astrophysics codes. Currently supported: X -> [Webnucleo](https://webnucleo.readthedocs.io) reaction XML and X -> [REACLIB](https://reaclib.jinaweb.org).

The idea: put your measured rate(T9) into X once, then convert that same file into whichever downstream format you need (Webnucleo and Reaclib now; PRISM, FRUITY planned) without re-entering the data. The Python package is imported as `unicorrn` (lowercase, per Python convention).

## Setup

Requires Python 3 (developed and tested with Python 3.12).

```bash
git clone https://github.com/chirag167/UniCoRRn.git
cd UniCoRRn
python3.12 -m venv venv
./venv/bin/python -m pip install -r requirements.txt
```

All commands below run scripts with `./venv/bin/python`, so you don't need to activate the venv. There's also nothing to `pip install` for UniCoRRn itself: the scripts load the `unicorrn` package straight from `src/`.

`pynucastro` is only needed for `xml_to_reaclib.py --validate`. The conversion itself doesn't require it.

The `venv/` folder is not committed (it's in `.gitignore`). If you move or rename the repo folder after creating it, the venv's scripts (e.g. `venv/bin/pip`) stop working because they store the venv's original path. Delete `venv/` and repeat the last two setup commands to rebuild it.

## Pipeline

```
                                                     --[convert_to_webnucleo.py]--> webnucleo XML
your_rate.csv --[csv_to_xml.py]--> X(TAMU) XML  <
                                                     --[xml_to_reaclib.py]------> REACLIB text
```

`build_example.py` is not part of this pipeline — it's a fixed demo/regression check for `src/unicorrn/reaction.py` and doesn't take real data.

## Step 1: the input CSV

One CSV file = **one reaction's rate table only**. Reactants, products, and all metadata (Q-value, source, evaluator, date, status) are *not* read from the CSV — they're passed as command-line flags when you run `csv_to_xml.py`. This keeps the CSV to just the numbers you actually measured.

**If you have multiple reactions, you need one CSV per reaction** (see "Multiple reactions" below for how they combine into a single X(TAMU) file).

Format: a header row naming the columns, then one row per T9 grid point.

| Column | Required? | Meaning |
|---|---|---|
| `t9` | required | temperature in units of 10^9 K |
| `value` or `rate` | required (either name works) | the reaction rate at that `t9` |
| `sigma_low` | optional | lower uncertainty bound on the rate |
| `sigma_high` | optional | upper uncertainty bound on the rate |

Notes:
- Column order doesn't matter; column *names* (case-insensitive) do.
- Blank lines and lines starting with `#` are ignored — use `#` for a comment/provenance line at the top of the file.
- Default delimiter is comma; pass `--delimiter` to `csv_to_xml.py` if your file uses something else (e.g. whitespace or tabs).
- If `sigma_low`/`sigma_high` are omitted entirely, or a cell is left blank, they come through as absent (not zero) — the converters treat that as "no uncertainty reported," not "uncertainty is zero."

Example (`data/c12_pg_n13_rate.csv`):

```
# FABRICATED numbers, for testing only -- not a real measured rate.
t9,value,sigma_low,sigma_high
0.010,1.2e-22,1.0e-22,1.5e-22
0.100,3.4e-11,3.0e-11,3.9e-11
0.500,2.1e-05,1.9e-05,2.3e-05
1.000,4.5e-03,4.1e-03,4.9e-03
2.000,1.8e-01,1.6e-01,2.0e-01
5.000,3.3e+00,3.0e+00,3.6e+00
```

A minimal file with no uncertainty columns is also valid (`data/o16_ag_ne20_rate.csv`):

```
t9,value
0.100,7.0e-18
0.500,2.4e-08
1.000,5.1e-04
2.000,3.2e-01
```

## Step 2: the mechanisms file

Every reaction also needs a `--mechanisms` JSON file: a non-empty list of the **physical component(s) that sum to the reaction's total rate**. This isn't optional bookkeeping — it's what lets `xml_to_reaclib.py` (see below) fit the rate the way REACLIB itself does, using physically-fixed coefficients rather than an unconstrained fit that can diverge outside your table.

Important: a rate is often the sum of *simultaneous* contributions (e.g. a non-resonant background plus one or more narrow resonances, all contributing at every T9) — not different mechanisms for different T9 ranges. List every component that's physically present, not "whichever one applies in a given range."

| Field | Applies to | Required? | Notes |
|---|---|---|---|
| `type` | all | required | `charge_induced_nonresonant` \| `neutron_induced_nonresonant` \| `narrow_resonant` |
| `l` | both non-resonant types | optional | minimum orbital angular momentum; defaults to `0` |
| `er` | `narrow_resonant` only | **required, no default** | resonance energy in MeV — an independently measured quantity, not derivable from the rate table |
| `s0` | `charge_induced_nonresonant` only | optional | astrophysical S-factor (MeV·barn) at zero energy; if omitted, `a0` is floated in the fit instead of fixed |

Example (one resonance plus a non-resonant background, both contributing):

```json
[{"type": "narrow_resonant", "er": 2.42},
 {"type": "charge_induced_nonresonant", "l": 1}]
```

A single-mechanism reaction is just a one-element list.

## Step 3: CSV -> X(TAMU)

```bash
./venv/bin/python csv_to_xml.py <rate_file.csv> <reaction_name> \
  --reactants <r1,r2> --products <p1,p2> \
  --mechanisms <mechanisms.json> \
  [--q <Q-value in MeV>] \
  [--source "<text>"] [--evaluator "<name>"] [--date <YYYY-MM-DD>] \
  [--status <recommended|upper_limit|placeholder|...>] \
  [--delimiter <char>] [--lab <lab tag, default TAMU>] \
  -o <output.xml> [--append] [--validate]
```

| Argument | Meaning |
|---|---|
| `rate_file` | the CSV described above |
| `name` | a unique item name for this reaction, e.g. `c12_pg_n13` — must be unique within the output file |
| `--reactants` / `--products` | comma-separated species names, e.g. `c12,p` / `n13,gamma` (your own naming — X(TAMU) doesn't enforce a convention) |
| `--mechanisms` | the mechanisms JSON file described above (required) |
| `--q` | Q-value in MeV (optional) |
| `--source`, `--evaluator`, `--date`, `--status` | free-text provenance metadata (all optional) |
| `-o, --output` | the X(TAMU) XML file to write |
| `--append` | add this reaction into an existing file at `-o` instead of overwriting it — this is how you combine reactions (see below) |
| `--validate` | check the result against xmlcoll's real schema (needs network access, since the schema imports a types file by URL) |

### Multiple reactions in one X(TAMU) file

Run the command once per reaction/CSV, reusing the same `-o` file with `--append` from the second call onward:

```bash
./venv/bin/python csv_to_xml.py data/c12_pg_n13_rate.csv c12_pg_n13 \
  --reactants c12,p --products n13,gamma \
  --mechanisms mechanisms_c12.json \
  --q 1.9435 --source "TAMU placeholder run" --evaluator "C. Rathi" \
  --date 2026-09-03 --status placeholder \
  -o data/cli_test.xml --validate

./venv/bin/python csv_to_xml.py data/o16_ag_ne20_rate.csv o16_ag_ne20 \
  --reactants o16,he4 --products ne20,gamma \
  --mechanisms mechanisms_o16.json \
  --status placeholder \
  -o data/cli_test.xml --append --validate
```

`data/cli_test.xml` now holds both reactions as separate items in one xmlcoll collection. Trying to `--append` a `name` that already exists in the file is an error (rename the reaction or fix the duplicate).

## Step 4: X(TAMU) -> Webnucleo

```bash
./venv/bin/python convert_to_webnucleo.py <x_tamu_file.xml> -o <output.xml> [--validate]
```

This reads *every* reaction item in the X(TAMU) file and writes them all into one webnucleo `reaction_data` XML file, using `wnutils`' own `Reaction`/`New_Xml` classes (no hand-built XML), then (with `--validate`) checks the result against libnucnet's real `reaction_data` schema offline.

```bash
./venv/bin/python convert_to_webnucleo.py data/cli_test.xml -o WebNucleo/webnucleo_reac.xml --validate
```

Two conversions happen automatically that are worth knowing about:

1. **Species names are translated** to webnucleo's nuclide-name convention (element symbol + mass number). Common shorthand is mapped automatically: `p`->`h1`, `d`->`h2`, `t`->`h3`, `a`/`alpha`->`he4`, `n`->`n1`, `g`->`gamma`, `e-`->`electron`, `e+`->`positron`. Anything else (e.g. `c12`, `ne20`) passes through unchanged — so use webnucleo-style names directly in `--reactants`/`--products` if you'd rather skip the alias table.
2. **Uncertainty is not carried over.** Webnucleo's native rate-table type only stores `(t9, rate, sef)` — there's no field for `sigma_low`/`sigma_high`. If your X(TAMU) reaction has uncertainties, the converter prints a warning and folds a note into the webnucleo file's `<source>` text, but the numbers themselves only remain machine-readable in the original X(TAMU) file. Keep that file as your record of uncertainty.
3. **The mechanisms field is silently unused here, on purpose.** Mechanism classification (`--mechanisms`) is fitting input specific to the REACLIB conversion (see below) — it isn't a universal fact about the rate the way uncertainty is, so webnucleo not using it doesn't lose anything and doesn't warrant a warning.

## Step 5: X(TAMU) -> Reaclib

```bash
./venv/bin/python xml_to_reaclib.py <x_tamu_file.xml> -o <output.reaclib> [--label <4-char>] [--validate]
```

```bash
./venv/bin/python xml_to_reaclib.py data/cli_test.xml -o Reaclib/reaclib_test.reaclib --validate
```

Unlike Webnucleo (which can hold your rate table directly), REACLIB's format only stores an analytic 7-coefficient exponential fit per reaction — there's no tabulated-rate option at all. So this conversion is a genuine **fit**, not a reshape, and it's done the way REACLIB itself does it: following the physically-constrained rules in Cyburt et al. (2010, the JINA REACLIB database paper, Table 1), most of the 7 coefficients (`a0..a6`) are *fixed* from the reaction's actual physics — charges, reduced mass, orbital angular momentum, resonance energy — and only the rest are floated against your table. That's exactly what `--mechanisms` (Step 2) supplies, and it's why this converter requires it while `convert_to_webnucleo.py` doesn't.

A naive fit that floats all 7 coefficients freely can fit a small table perfectly (0% residual) yet diverge to infinity just outside the tabulated range, because nothing constrains its asymptotic behavior — fixing the physically-determined coefficients is what prevents that.

A few behaviors worth knowing about:

1. **Every reaction needs `mechanisms` and a Q-value.** Both are required by REACLIB's format and can't be defaulted — a reaction missing either is skipped with a clear error (the rest of the file still converts).
2. **Only 2-reactant reactions are supported** (REACLIB chapters 4, 5, 6, 7). Table 1's fitting rules describe two-body entrance-channel physics; 1- or 3-reactant reactions (e.g. decays, triple-alpha) aren't implemented yet and are skipped with an error rather than silently mishandled.
3. **Emitted photons are dropped from the product list, not translated.** REACLIB doesn't track a photon as a nuclide token at all — `c12(p,g)n13` is chapter 4 with tokens `[c12, p, n13]`, not chapter 5 with an explicit `gamma`/`g` product. The converter strips `gamma`/`g`/`photon` from products before computing the chapter and writing the file.
4. **Fits exceeding REACLIB's own 5% accuracy target are flagged, not hidden or rejected.** REACLIB's own paper states a proper fit should be within 5% of the data; this converter checks the same threshold and prints a `DEVIATION` warning to stderr for any reaction that misses it — mirroring REACLIB's own "deviations list" practice — but still writes it, since 5% may simply be unreachable with sparse data.
5. **Multiple mechanisms are summed and fit jointly**, which is a genuinely nonlinear fit (a sum of exponentials isn't linear in the coefficients the way a single set is), using `scipy.optimize.least_squares`.
6. **A 0% residual is not automatically good news.** If a reaction has fewer tabulated points than floated parameters, the fit is underdetermined and will interpolate exactly through every point regardless of whether it represents the data at all — the same trap the very first unconstrained-fit prototype hit. Check the point-count-vs-floated-parameter balance before trusting a suspiciously perfect fit.
7. **`--validate`** reads the output back with pynucastro's real `ReacLibLibrary` reader and cross-checks its evaluated rate against this converter's own in-memory fit — confirming the ASCII serialization round-trips correctly (not re-checking fit quality, which is reported separately per point 4).

## Output conventions in this repo

- `data/` — input CSVs, mechanisms JSON files, and X(TAMU) XML files
- `WebNucleo/` — converted webnucleo XML files (point `convert_to_webnucleo.py -o` here directly)
- `Reaclib/` — converted REACLIB text files (point `xml_to_reaclib.py -o` here directly)

## Files

| File | Role |
|---|---|
| `src/unicorrn/reaction.py` | shared read/write logic (`build_reaction_item`, `read_reaction_item`) all CLI scripts build on |
| `src/unicorrn/reaclib_fit.py` | REACLIB chapter lookup, Table 1's fixed/float coefficient rules, and the joint nonlinear fit |
| `csv_to_xml.py` | CLI: CSV -> X(TAMU) XML |
| `convert_to_webnucleo.py` | CLI: X(TAMU) XML -> webnucleo `reaction_data` XML |
| `xml_to_reaclib.py` | CLI: X(TAMU) XML -> REACLIB text |
| `build_example.py` | fixed demo/regression check, not part of the real pipeline |
