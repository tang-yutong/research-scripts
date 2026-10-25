#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""phd-radar: a small daily radar for PhD vacancies.

Scans a few official, publicly available job boards (ETH Zurich jobs portal,
AcademicTransfer for the Netherlands, Euraxess for the EU), keeps a local
JSON state file so every run reports only what is NEW since the last run,
filters postings with keyword rules, and writes a Markdown digest.

The parsing logic is factored into pure functions (parse_ethz / parse_at /
parse_euraxess / should_keep / merge_seen) so it can be unit-tested against
synthetic HTML without any network access.

Etiquette: one page per source/keyword per run, plain stdlib requests with a
browser UA -- a daily personal digest, not a crawler. Standard library only.
"""
import argparse
import datetime
import html
import json
import pathlib
import re
import sys
import urllib.parse
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from config import CONFIG  # noqa: E402
import common  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")

UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                    "(KHTML, like Gecko) Chrome/129.0 Safari/537.36"}

# ================== filtering rules -- edit to taste (examples) ==================
# a posting title must look like a doctoral-level opening
PHD = re.compile(r"PhD|doctoral|doktorand|predoc|doctorant|R1", re.I)
# keep if the title matches any keyword (example: a chemistry-oriented set)
RELEVANT = re.compile(
    r"chem|molecular|molecul|simulat|DFT|quantum|comput|materials?|catalys|"
    r"theoret|dynamics|drug|protein|biophys", re.I)
# drop if the title matches any of these (example starter set -- extend freely)
EXCLUDE = re.compile(r"nursing|architecture|phonetics|theology", re.I)
# per-board search terms (edit to your field)
# dedup state keeps at most this many seen keys (oldest dropped first);
# the authoritative value lives in config.radar_state_cap
AT_QUERIES = ("chemistry", "computational", "molecular simulation", "catalysis")
EURAXESS_QUERIES = ("computational chemistry", "molecular dynamics", "DFT",
                    "catalysis", "medicinal chemistry", "drug discovery")
# =================================================================================


def should_keep(title: str, phd: re.Pattern = PHD, relevant: re.Pattern = RELEVANT,
                exclude: re.Pattern = EXCLUDE) -> bool:
    """True iff a title passes all three keyword filters."""
    return bool(phd.search(title) and relevant.search(title)
                and not exclude.search(title))


def merge_seen(prev, new_keys, cap=None):
    """Append new keys to the seen list, drop duplicates, keep the newest
    `cap` entries (sliding window; default from config.radar_state_cap)."""
    if cap is None:
        cap = CONFIG.radar_state_cap
    combined = list(dict.fromkeys(list(prev) + list(new_keys)))
    return combined[-cap:]


def parse_ethz(html_text: str):
    """ETH jobs portal (server-rendered list) -> [posting dicts]."""
    got = {}
    for jid, inner in re.findall(
            r'<a[^>]+href="/job/view/(JOPG_ethz_[A-Za-z0-9_]+)"[^>]*>(.*?)</a>',
            html_text, re.S):
        title = re.sub(r"\s+", " ",
                       html.unescape(re.sub(r"<[^>]+>", " ", inner))).strip()
        got["ethz:" + jid] = {"dl": "?", "title": title[:90],
                              "src": "jobs.ethz.ch/job/view/" + jid,
                              "key": "ethz:" + jid, "org": "ETH Zurich"}
    return list(got.values())


def parse_at(html_text: str):
    """AcademicTransfer search page -> [posting dicts]."""
    got = {}
    for _href, jid, slug in re.findall(r'href="(/en/jobs/(\d+)/([^"]+))"',
                                       html_text):
        title = slug.replace("-", " ").strip()
        got["at:" + jid] = {"dl": "?", "title": title[:90],
                            "src": "academictransfer.com/en/jobs/" + jid,
                            "key": "at:" + jid, "org": "AcademicTransfer"}
    return list(got.values())


def parse_euraxess(html_text: str):
    """Euraxess keyword-search page -> [posting dicts]; cards are anchored by
    /jobs/<id>, title = first span after the anchor, deadline = first <time>."""
    out = []
    anchors = [(m.start(), m.group(1)) for m in re.finditer(r'href="/jobs/(\d+)"',
                                                            html_text)]
    done = set()
    for i, (pos, jid) in enumerate(anchors):
        if jid in done:
            continue
        done.add(jid)
        seg = html_text[pos:anchors[i + 1][0] if i + 1 < len(anchors)
                        else pos + 8000]
        tm = re.search(r'href="/jobs/' + jid +
                       r'"[^>]*>\s*(?:<span>)?\s*([^<]{8,200})', seg)
        title = re.sub(r"\s+", " ", html.unescape(tm.group(1))).strip() if tm else ""
        if "<" in title or "class=" in title or not title:
            title = "(title not parsed, open detail page) euraxess:" + jid
        dlm = re.search(r"<time[^>]*>([^<]+)</time>", seg)
        dl = re.sub(r"\s+", " ", dlm.group(1)).strip()[:22] if dlm else "?"
        out.append({"dl": dl, "title": title[:90], "src": "euraxess:" + jid,
                    "key": "eu:" + jid, "org": "Euraxess"})
    return out


def fetch(url: str, timeout: int = 40) -> str:
    req = urllib.request.Request(url, headers=UA)
    return urllib.request.urlopen(req, timeout=timeout).read().decode(
        "utf-8", "ignore")


def load_state(state_path: pathlib.Path) -> dict:
    if state_path.exists():
        return json.loads(state_path.read_text(encoding="utf-8"))
    return {"seen": []}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Daily PhD-vacancy radar.")
    common.add_common_args(ap)
    args = ap.parse_args(argv)
    log = common.get_logger("phd_radar", args.log_level)
    cfg = CONFIG
    if args.config:
        common.apply_overrides(cfg, args.config)
    common.resolve(cfg)

    state_path = cfg.out_dir / "radar_state.json"
    digest_path = cfg.out_dir / ("digest_" + datetime.date.today().isoformat()
                                 + ".md")
    if args.dry_run:
        log.info("would fetch %d queries across 3 boards; state %s; digest %s",
                 1 + len(AT_QUERIES) + len(EURAXESS_QUERIES), state_path,
                 digest_path)
        return common.EXIT_OK

    state = load_state(state_path)
    seen = list(state["seen"])
    known = set(seen)

    sources = []
    # each source is isolated: one unreachable board never kills the run
    try:
        sources += [p for p in parse_ethz(fetch("https://www.jobs.ethz.ch/"))
                    if should_keep(p["title"])]
    except Exception as e:
        sources.append({"dl": "-", "title": "[source unreachable] ETH portal: "
                                            + str(e)[:60], "src": "-",
                        "key": "ethz-err", "org": ""})
    errs = []
    for q in AT_QUERIES:
        try:
            sources += [p for p in parse_at(fetch(
                "https://www.academictransfer.com/en/jobs/?q="
                + urllib.parse.quote(q) + "&category=phd")) if should_keep(p["title"])]
        except Exception as e:
            errs.append(str(e)[:50])
    if errs and not any(p["key"].startswith("at:") for p in sources):
        sources.append({"dl": "-", "title": "[source unreachable] "
                        "AcademicTransfer: " + "; ".join(errs), "src": "-",
                        "key": "at-err", "org": ""})
    eu_errs = []
    for kw in EURAXESS_QUERIES:
        try:
            sources += [p for p in parse_euraxess(fetch(
                "https://euraxess.ec.europa.eu/jobs/search?f%5B0%5D=keywords%3A"
                + urllib.parse.quote(kw))) if should_keep(p["title"])]
        except Exception as e:
            eu_errs.append((kw, str(e)[:50]))
    for kw, err in eu_errs:
        sources.append({"dl": "-", "title": f"[source unreachable] Euraxess({kw}): {err}",
                        "src": "-", "key": "eu-err:" + kw, "org": ""})

    fresh = [x for x in sources if x["key"].endswith("-err")
             or x["key"].startswith("eu-err:") or x["key"] not in known]
    seen = merge_seen(seen, [x["key"] for x in sources if not x["key"].endswith(
        ("-err",)) and not x["key"].startswith("eu-err:")])
    state_path.parent.mkdir(parents=True, exist_ok=True)
    state_path.write_text(json.dumps({"seen": seen}, ensure_ascii=False),
                          encoding="utf-8")

    n_new = len([x for x in fresh if not x["key"].endswith("-err")
                 and not x["key"].startswith("eu-err:")])
    lines = ["# job digest " + datetime.date.today().isoformat(), "",
             f"> new this run: {n_new} (auto pre-filter: PhD-level + keyword "
             "rules; review manually before acting)", ""]
    if not fresh:
        lines.append("(nothing new)")
    for x in sorted(fresh, key=lambda r: r["dl"]):
        lines.append(f"- **{x['dl']}** | {x['title']} | {x['org']} | {x['src']}")
    digest_path.parent.mkdir(parents=True, exist_ok=True)
    digest_path.write_text("\n".join(lines), encoding="utf-8")
    log.info("wrote %s | new: %d", digest_path, n_new)
    return 0


if __name__ == "__main__":
    sys.exit(main())
