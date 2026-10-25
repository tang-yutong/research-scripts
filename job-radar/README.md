# job-radar

`phd_radar.py` -- a small daily radar for PhD vacancies. Scans three official,
publicly available boards (ETH Zurich jobs portal, AcademicTransfer for the
Netherlands, Euraxess for the EU), deduplicates against a local JSON state
file so every run reports only **new** postings, filters with keyword rules,
and writes a Markdown digest.

Standard library only: `python phd_radar.py`.

## Configure (top of the file, marked as examples)

- `PHD` / `RELEVANT` / `EXCLUDE` -- keyword filters; the shipped set is one
  field's example, edit to yours.
- `AT_QUERIES` / `EURAXESS_QUERIES` -- per-board search terms.

Paths (digest and state files land under `out/`, like every other module)
and the dedup window size (`radar_state_cap`) live in `config.py`; the
state is a sliding window -- oldest seen keys drop out beyond the cap.
Both files are local only and gitignored via `out/`.

## Testability

The board parsers and the filters are pure functions -- `parse_ethz`,
`parse_at`, `parse_euraxess`, `should_keep`, `merge_seen` -- covered by
`tests/test_radar.py` against synthetic HTML, no network needed.

## Etiquette

One page per source/keyword per run, plain stdlib requests with a browser
UA -- a daily personal digest, not a crawler. Keep it that way.
