"""Tests for job-radar: keyword filters, dedup-state merge, board parsers
(all against synthetic HTML -- no network)."""
from phd_radar import should_keep, merge_seen, parse_ethz, parse_at, parse_euraxess

ETHZ_HTML = """
<div>
<a href="/job/view/JOPG_ethz_ABC123">Doctoral Student in Computational Chemistry 100%</a>
<p>some details</p>
<a href="/job/view/JOPG_ethz_DEF456">Senior Lecturer in Musicology</a>
</div>
"""

AT_HTML = """
<a href="/en/jobs/123456/doctoral-candidate-in-molecular-simulation">x</a>
<a href="/en/jobs/123457/professor-of-architectural-history">y</a>
"""

EURAXESS_HTML = """
<a href="/jobs/472218"><span>PhD Candidate in Molecular Dynamics Simulations</span></a>
<div><time>1 Nov 2026 - 22:59 (UTC)</time> Brussels, Belgium</div>
<a href="/jobs/474960"><span>Doctoral Researcher in solid-state NMR</span></a>
<div><time>15 Nov 2026 - 22:59 (UTC)</time></div>
"""


def test_radar_should_keep():
    assert should_keep("Doctoral Candidate in Computational Chemistry")
    assert should_keep("PhD position in molecular dynamics")
    assert not should_keep("Postdoc in NMR spectroscopy")        # not PhD-level
    assert not should_keep("PhD in Architectural History")       # excluded field
    assert not should_keep("Doctoral position in medieval art")  # not relevant


def test_radar_merge_seen_cap():
    seen = merge_seen(["a", "b"], ["b", "c"])
    assert seen == ["a", "b", "c"]                # deduped, order kept
    assert merge_seen([], ["x", "y"]) == ["x", "y"]
    out = merge_seen(["k%d" % i for i in range(10)], ["new"], cap=5)
    assert out == ["k6", "k7", "k8", "k9", "new"]  # oldest dropped first


def test_radar_parse_ethz():
    postings = parse_ethz(ETHZ_HTML)
    assert len(postings) == 2
    by_key = {p["key"]: p for p in postings}
    assert "ethz:JOPG_ethz_ABC123" in by_key
    assert "Doctoral Student in Computational Chemistry" in by_key[
        "ethz:JOPG_ethz_ABC123"]["title"]
    # parsing is filter-agnostic; the filter separates them
    keep = [p for p in postings if should_keep(p["title"])]
    assert [p["key"] for p in keep] == ["ethz:JOPG_ethz_ABC123"]


def test_radar_parse_at():
    postings = parse_at(AT_HTML)
    by_key = {p["key"]: p for p in postings}
    assert by_key["at:123456"]["title"] == "doctoral candidate in molecular simulation"
    keep = [p for p in postings if should_keep(p["title"])]
    assert [p["key"] for p in keep] == ["at:123456"]


def test_radar_parse_euraxess():
    postings = parse_euraxess(EURAXESS_HTML)
    assert len(postings) == 2
    by_key = {p["key"]: p for p in postings}
    assert by_key["eu:472218"]["dl"].startswith("1 Nov 2026")
    assert "Molecular Dynamics" in by_key["eu:472218"]["title"]
    assert by_key["eu:474960"]["org"] == "Euraxess"
