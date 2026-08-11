#!/usr/bin/env python3
"""
check_dois.py -- validate every DOI in EmergingXOR.tex against Crossref.

Twenty-six references were added during revision and written from memory.  A
spot-check of seven found one wrong (a journal issue number), so the expected
error rate in the remainder is not negligible and a full check is worth the two
minutes it takes.  This script extracts each DOI from the .tex, asks Crossref
what that DOI actually is, and flags any entry whose year, volume, issue, pages
or first author disagree with what the manuscript claims.

It reports rather than edits: bibliographic corrections should be made by hand.

    python3 check_dois.py                     # defaults to EmergingXOR.tex
    python3 check_dois.py --file EmergingXOR_rev.tex
    python3 check_dois.py --verbose           # also print matches

Needs network access to api.crossref.org.  Put a real address in MAILTO: the
Crossref "polite pool" is faster and it is common courtesy.
"""
import argparse, json, re, sys, time, unicodedata
import urllib.parse
import urllib.request

VERSION = "1.2"   # 1.0: false mismatches from LaTeX accents, Unicode
                  # hyphens, name particles, unescaped \_ in DOIs.
                  # 1.1: page regex swallowed the sentence-final period, and
                  # the loop variable shadowed the consistency counter.
MAILTO = "marco@ifca.unican.es"
API = "https://api.crossref.org/works/"

DOI_RE = re.compile(r"https?://doi\.org/(10\.[^\s}]+)")
# \item Author, A. B., & Other, C. (2016). Title. \emph{Journal, 36}(40), 1--2. ...
# \Z is required: the slice below stops before \end{apabib}, so without it the
# final reference is silently dropped.
ITEM_RE = re.compile(r"\\item\s+(.+?)(?=\n\\item|\n\\end\{apabib\}|\Z)", re.S)


ACCENTS = {r"\\'": "", r'\\"': "", r"\\`": "", r"\\\^": "", r"\\~": "",
           r"\\c": "", r"\\v": "", r"\\=": "", r"\\\.": ""}


def norm(s):
    s = re.sub(r"\\[a-zA-Z]+\{([^}]*)\}", r"\1", s)      # strip \emph{...}
    # LaTeX accents: \"o -> o, \'e -> e, \~n -> n, also the braced forms \'{e}
    s = re.sub(r"\\[\'`\"^~=.]\{?([a-zA-Z])\}?", r"\1", s)
    s = re.sub(r"\\[a-zA-Z]", "", s)
    s = re.sub(r"[\\{}$~]", "", s)
    s = re.sub(r"\s+", " ", s)
    return s.strip()


def fold(s):
    """Compare names accent- and hyphen-insensitively.

    Crossref returns composed Unicode ('Bohm' with an umlaut) and sometimes the
    non-breaking hyphen U+2010, while the .tex carries LaTeX escapes and ASCII
    hyphens.  Neither difference is a bibliographic error, so both sides are
    folded to bare ASCII before comparison.
    """
    s = unicodedata.normalize("NFKD", str(s))
    s = "".join(c for c in s if not unicodedata.combining(c))
    for dash in "\u2010\u2011\u2012\u2013\u2014":
        s = s.replace(dash, "-")
    return s.lower().strip()


def same_name(mine, theirs):
    """Surnames match if either contains the other once folded: the manuscript
    may drop a particle ('De Angeli' -> 'Angeli') without being wrong."""
    a, b = fold(mine), fold(theirs)
    return a == b or a in b or b in a


def fetch(doi):
    url = API + urllib.parse.quote(doi) + f"?mailto={MAILTO}"
    req = urllib.request.Request(url, headers={"User-Agent": f"check_dois/1.0 (mailto:{MAILTO})"})
    with urllib.request.urlopen(req, timeout=20) as fh:
        return json.load(fh)["message"]


def trim_pages(first, last):
    """Join a page range, dropping the sentence-final period.

    The page pattern has to allow '.' so that elocators like '655-664.e4'
    survive, which means it also swallows the full stop that ends the entry.
    Strip trailing dots but keep an internal one.
    """
    return f"{first.strip('.')}-{last.rstrip('.')}"


def claimed(entry):
    """Pull year / volume / issue / pages as the manuscript states them."""
    out = {}
    m = re.search(r"\((\d{4})[a-c]?\)", entry)
    if m:
        out["year"] = m.group(1)
    # \emph{Journal, VOL}(ISS), PAGES
    m = re.search(r",\s*(\d+)\}\((\d+)\),\s*([\dA-Za-z.]+)(?:--|–)([\dA-Za-z.]+)", entry)
    if m:
        out["volume"], out["issue"] = m.group(1), m.group(2)
        out["pages"] = trim_pages(m.group(3), m.group(4))
    else:
        m = re.search(r",\s*(\d+)\},\s*([\dA-Za-z.]+)(?:--|–)([\dA-Za-z.]+)", entry)
        if m:
            out["volume"] = m.group(1)
            out["pages"] = trim_pages(m.group(2), m.group(3))
    m = re.match(r"([A-ZÀ-Ÿ][^,]*),", norm(entry))
    if m:
        out["surname"] = m.group(1).split()[-1]
    return out


def actual(msg):
    out = {}
    for key in ("published-print", "published-online", "issued"):
        parts = msg.get(key, {}).get("date-parts", [[None]])[0]
        if parts and parts[0]:
            out["year"] = str(parts[0]); break
    for a, b in (("volume", "volume"), ("issue", "issue"), ("page", "pages")):
        if msg.get(a):
            out[b] = str(msg[a]).replace("–", "-")
    auth = msg.get("author") or []
    if auth and auth[0].get("family"):
        out["surname"] = auth[0]["family"]
    out["title"] = (msg.get("title") or [""])[0]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--file", default="EmergingXOR.tex")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--sleep", type=float, default=0.3)
    a = ap.parse_args()

    print(f"check_dois.py v{VERSION}")
    tex = open(a.file, encoding="utf-8").read()
    try:
        block = tex[tex.index(r"\begin{apabib}"):tex.index(r"\end{apabib}")]
    except ValueError:
        print("no apabib block found"); sys.exit(1)

    entries = ITEM_RE.findall(block)
    print(f"{a.file}: {len(entries)} reference entries\n")

    bad = miss = ok = 0
    for e in entries:
        m = DOI_RE.search(e)
        label = norm(e)[:60]
        if not m:
            if "arxiv" in e.lower() or "Patent" in e:
                continue                      # arXiv preprints and patents: no DOI
            print(f"[NO DOI ] {label}"); miss += 1; continue
        doi = m.group(1).rstrip(".,;").replace("\\_", "_").replace("\\&", "&")
        try:
            msg = fetch(doi)
        except Exception as exc:
            print(f"[FAIL   ] {label}\n           {doi} -> {exc}"); miss += 1
            time.sleep(a.sleep); continue
        c, r = claimed(e), actual(msg)
        diffs = []
        for k in ("year", "volume", "issue", "pages", "surname"):
            if k not in c or not r.get(k):
                continue
            agree = same_name(c[k], r[k]) if k == "surname" else fold(c[k]) == fold(r[k])
            if not agree:
                diffs.append((k, c[k], r[k]))
        if diffs:
            print(f"[MISMATCH] {label}")
            print(f"           doi   : {doi}")
            print(f"           actual: {r.get('title','')[:70]}")
            for k, mine, theirs in diffs:
                print(f"           {k:8s} manuscript={mine!r}  crossref={theirs!r}")
            bad += 1
        else:
            ok += 1
            if a.verbose:
                print(f"[OK      ] {label}")
        time.sleep(a.sleep)

    print(f"\n{ok} consistent, {bad} mismatched, {miss} unresolved")
    if bad:
        print("Mismatches are usually a wrong issue or page range, not a wrong paper;")
        print("check each against the DOI landing page before editing.")


if __name__ == "__main__":
    main()
