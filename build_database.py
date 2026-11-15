
"""Tissue Regeneration Peptide Navigator - Europe PMC builder v2.

Original Excel and existing SQLite records are preserved.
Run:
 python build_database.py --self-test
 python build_database.py --report
 python build_database.py --peptide GHK --limit 10
 python build_database.py --all --limit 10
 python build_database.py --peptide GHK --limit 10 --rescreen-existing
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
import time
import xml.etree.ElementTree as ET

from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import pandas as pd

ROOT = Path(__file__).resolve().parent
EXCEL = ROOT / "peptide_regeneration_DB_3차_screening.xlsx"
DB = ROOT / "navigator.sqlite3"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest"
NA = "Not reported"

ALIASES = {
    "GHK": ['"GHK"', '"glycyl-histidyl-lysine"', '"GHK-Cu"'],
    "RADA16-I": ['"RADA16"', '"RADA16-I"'],
    "QK": ['"QK peptide"', '"VEGF mimetic peptide"'],
    "IKVAV": ['"IKVAV"'],
    "KLD-12": ['"KLD12"', '"KLD-12"'],
    "RGD": ['"RGD peptide"', '"arginine-glycine-aspartate"'],
    "YIGSR": ['"YIGSR"'],
}

MENTIONS = {
    "GHK": r"(?<![A-Za-z])(?:GHK(?:[- ](?:Cu|C16))?|Cu[- ]GHK(?:[- ]C16)?|glycyl[- ]histidyl[- ]lysine)(?![A-Za-z])",
    "RADA16-I": r"(?<![A-Za-z])RADA16(?:[- ](?:I|RGD|IKVAV))?(?![A-Za-z])",
    "QK": r"\b(?:QK\s+peptide|QK[- ](?:GelMA|DA)|VEGF mimetic peptide)\b",
    "IKVAV": r"\bIKVAV\b",
    "KLD-12": r"\bKLD[- ]?12\b",
    "RGD": r"\b(?:RGD\s+(?:peptide|motif|sequence)|RADA16[- ]RGD)\b",
    "YIGSR": r"\bYIGSR\b",
}

FUNCTIONS = {
    "WOUND": r"wound|epitheli|keratinocyte|fibroblast|scratch|healing",
    "ANGIO": r"angiogen|vascular|endothelial|tube formation|CD31",
    "NERVE": r"neur|axon|schwann|myelin|sciatic|neural",
    "CARTILAGE": r"cartilage|chondrogen",
    "ADHESION": r"cell adhesion|cell attachment|cell spreading",
    "ECMREMOD": r"collagen|extracellular matrix|ECM remodel",
}

OUTCOME = re.compile(
    r"\b(?:increas\w*|decreas\w*|enhanc\w*|reduc\w*|"
    r"promot\w*|improv\w*|inhibit\w*|significan\w*|"
    r"compar\w*|higher|lower|migration|regenerat\w*)\b",
    re.I,
)

REVIEW = re.compile(
    r"\b(?:review|meta[- ]analysis|bibliometric|"
    r"systematic review|scoping review)\b",
    re.I,
)

CITATION = re.compile(
    r"\[\s*\d+(?:\s*[,–-]\s*\d+)*\s*\]|"
    r"\b[A-Z][a-z]+\s+et al\.?\s*\[?\d+\]?"
)

SENT_SPLIT = re.compile(r"(?<=[.!?])\s+(?=[A-Z(])")


def clean(x):
    return re.sub(r"\s+", " ", str(x or "")).strip()


def now():
    return datetime.now(timezone.utc).isoformat()


def doi_norm(x):
    return re.sub(
        r"^https?://(?:dx\.)?doi\.org/",
        "",
        clean(x),
        flags=re.I,
    ).lower()


def get(url, params=None, as_xml=False):
    if params:
        url += "?" + urlencode(params)

    request = Request(
        url,
        headers={
            "User-Agent": "TissuePeptideNavigator/2.0 (research)",
            "Accept": (
                "application/xml" if as_xml
                else "application/json"
            ),
        },
    )

    for attempt in range(3):
        try:
            with urlopen(request, timeout=25) as resp:
                raw = resp.read(12_000_001)

            if len(raw) > 12_000_000:
                raise ValueError("API response exceeds 12MB")

            return (
                ET.fromstring(raw)
                if as_xml else json.loads(raw)
            )

        except HTTPError as exc:
            if exc.code in (400, 401, 403, 404) or attempt == 2:
                raise

        except (URLError, TimeoutError):
            if attempt == 2:
                raise

        time.sleep(2 ** attempt)


def initialize(db):
    db.executescript("""
    PRAGMA foreign_keys=ON;

    CREATE TABLE IF NOT EXISTS excel_screening (
        record_id TEXT PRIMARY KEY,
        peptide TEXT,
        function_code TEXT,
        study_id TEXT,
        doi TEXT,
        payload_json TEXT NOT NULL
    );

    CREATE TABLE IF NOT EXISTS excel_recommendations (
        recommendation_id TEXT PRIMARY KEY,
        peptide_id TEXT,
        function_code TEXT,
        doi TEXT,
        payload_json TEXT NOT NULL
    );

    CREATE TABLE IF NOT EXISTS papers (
        paper_key TEXT PRIMARY KEY,
        title TEXT,
        year TEXT,
        doi TEXT,
        pmid TEXT,
        pmcid TEXT,
        paper_url TEXT,
        abstract TEXT,
        fulltext_status TEXT,
        collected_at TEXT
    );

    CREATE TABLE IF NOT EXISTS paper_peptides (
        paper_key TEXT NOT NULL REFERENCES papers(paper_key),
        peptide TEXT NOT NULL,
        PRIMARY KEY(paper_key, peptide)
    );

    CREATE TABLE IF NOT EXISTS extraction_candidates (
        candidate_id TEXT PRIMARY KEY,
        paper_key TEXT NOT NULL REFERENCES papers(paper_key),
        peptide TEXT,
        source_section TEXT,
        source_text TEXT NOT NULL,
        screening_json TEXT NOT NULL,
        evidence_status TEXT NOT NULL,
        created_at TEXT NOT NULL
    );

    CREATE INDEX IF NOT EXISTS idx_candidate_peptide
    ON extraction_candidates(peptide);

    CREATE TABLE IF NOT EXISTS paper_screening_v2 (
        paper_key TEXT PRIMARY KEY,
        study_type TEXT NOT NULL,
        classification_reason TEXT,
        metadata_json TEXT,
        assessment_at TEXT NOT NULL
    );

    CREATE TABLE IF NOT EXISTS candidate_context_v2 (
        candidate_id TEXT PRIMARY KEY,
        paper_key TEXT NOT NULL,
        peptide TEXT NOT NULL,
        section_title TEXT,
        paragraph TEXT,
        prior_paragraph TEXT,
        mention_scope TEXT,
        citation_warning INTEGER,
        is_primary_candidate INTEGER NOT NULL DEFAULT 0,
        collected_at TEXT NOT NULL
    );
    """)

    db.commit()


def import_excel(db, path):
    sheets = pd.read_excel(
        path,
        sheet_name=[
            "Screening DB v3",
            "Recommendation DB v3",
            "Peptide Library",
        ],
        dtype=str,
        engine="openpyxl",
        keep_default_na=False,
    )

    available = {
        clean(x) for x in sheets["Peptide Library"]["Name"]
    }

    if set(ALIASES) - available:
        raise ValueError(
            "Peptide Library missing: "
            + str(sorted(set(ALIASES) - available))
        )

    screening = sheets["Screening DB v3"]
    recommendations = sheets["Recommendation DB v3"]

    with db:
        for _, row in screening.iterrows():
            item = row.to_dict()

            db.execute(
                "INSERT OR IGNORE INTO excel_screening "
                "VALUES (?,?,?,?,?,?)",
                (
                    clean(item["Record_ID"]),
                    clean(item["Original_Peptide"]),
                    clean(item["Target_Function_Code"]),
                    clean(item["Study_ID"]),
                    doi_norm(item["DOI"]),
                    json.dumps(item, ensure_ascii=False),
                ),
            )

        for _, row in recommendations.iterrows():
            item = row.to_dict()

            db.execute(
                "INSERT OR IGNORE INTO excel_recommendations "
                "VALUES (?,?,?,?,?)",
                (
                    clean(item["Recommendation_ID"]),
                    clean(item["Peptide_ID"]),
                    clean(item["Target_Function_Code"]),
                    doi_norm(item["DOI"]),
                    json.dumps(item, ensure_ascii=False),
                ),
            )

    return len(screening), len(recommendations)


def identity(p):
    if doi_norm(p.get("doi")):
        return "DOI:" + doi_norm(p["doi"])

    if clean(p.get("pmid")):
        return "PMID:" + clean(p["pmid"])

    if clean(p.get("pmcid")):
        return "PMCID:" + clean(p["pmcid"])

    return ""


def classify_publication(p, xml_root=None):
    title = clean(p.get("title"))

    types = p.get("pubTypeList") or {}

    if isinstance(types, dict):
        types = types.get("pubType", [])

    if isinstance(types, str):
        types = [types]

    types = [clean(t) for t in types]

    if any(REVIEW.search(t) for t in types):
        return (
            "review",
            "Europe PMC publication type: " + ", ".join(types),
        )

    if REVIEW.search(title):
        return "review_possible", "Review indicator in title"

    if xml_root is not None:
        kind = xml_root.attrib.get("article-type", "").lower()

        if "review" in kind:
            return "review", "JATS article-type=" + kind

        has_results = bool(
            xml_root.findall(".//body//sec")
        )

        if (
            has_results
            and kind in ("research-article", "original-article")
        ):
            return (
                "original_candidate",
                "JATS article-type=" + kind,
            )

    return (
        "unknown",
        "Article type not confirmed as original research",
    )


def query_for(peptide):
    aliases = " OR ".join(ALIASES[peptide])

    return (
        "(" + aliases + ") "
        "AND (wound OR regenerat* OR angiogen* "
        "OR nerve OR tissue OR hydrogel) "
        "AND (SRC:MED OR SRC:PMC)"
    )


def content_blocks(root):
    body = root.find(".//body")

    if body is None:
        return []

    blocks = []

    def walk(sec, inherited=""):
        heading = sec.find("title")

        name = (
            clean(" ".join(heading.itertext()))
            if heading is not None else inherited
        )

        tag = (
            (inherited + " / " + name).strip(" /")
            if name != inherited else inherited
        )

        is_results = bool(
            re.search(
                r"\bresults?\b|\bfindings?\b|results and discussion",
                tag,
                re.I,
            )
        )

        previous = ""

        if is_results:
            for para in sec.findall("./p"):
                text = clean(" ".join(para.itertext()))

                if text:
                    blocks.append((tag or "Results", text, previous))
                    previous = text

        for child in sec.findall("./sec"):
            walk(child, tag)

    for section in body.findall("./sec"):
        walk(section)

    return blocks


def make_candidates(paragraph, previous, peptide, section):
    results = []

    for sentence in SENT_SPLIT.split(paragraph):
        sentence = clean(sentence)

        if not 35 <= len(sentence) <= 1200:
            continue

        if not OUTCOME.search(sentence):
            continue

        mention = bool(
            re.search(MENTIONS[peptide], sentence, re.I)
        )

        if not mention:
            continue

        matched = [
            code
            for code, pattern in FUNCTIONS.items()
            if re.search(pattern, sentence, re.I)
        ]

        if not matched:
            continue

        citation_detected = bool(
            CITATION.search(sentence)
        )

        for code in matched:
            payload = {
                "Record_ID": NA,
                "Original_Peptide": peptide,
                "Original_Sequence": NA,
                "Modification_Position": NA,
                "Modification_Method": NA,
                "Modified_Peptide": NA,
                "Modified_Sequence": NA,
                "Formulation": NA,
                "Intervention_Class": "확인 필요",
                "Target_Function": NA,
                "Target_Function_Code": code,
                "Evaluation_Method": NA,
                "Measured_Outcome": NA,
                "Effect_Direction": "확인 필요",
                "Quantitative_Result": NA,
                "Statistical_Significance": NA,
                "Comparator": NA,
                "Comparator_Type": NA,
                "Evidence_Status": "UNVERIFIED_EXTRACTION",
                "Source_Verification": section,
                "Reference": NA,
                "DOI": NA,
                "Paper_URL": NA,
                "Study_ID": NA,
                "Candidate_Statement": sentence,
                "Reviewer_Notes": (
                    "Results candidate only; "
                    "context preserved separately. "
                    "Manual or stricter verification required."
                ),
            }

            results.append(
                (sentence, code, payload, int(citation_detected))
            )

    return results


def save_paper(
    db,
    p,
    peptide,
    fulltext=True,
    rescreen_existing=False,
):
    key = identity(p)

    if not key:
        return 0, 0, "no_identifier"

    pmcid = clean(p.get("pmcid"))
    doi = doi_norm(p.get("doi"))
    abstract = clean(p.get("abstractText"))

    url = (
        "https://europepmc.org/articles/" + pmcid
        if pmcid else
        "https://doi.org/" + doi
        if doi else
        "https://europepmc.org/article/MED/"
        + clean(p.get("pmid"))
    )

    present = db.execute(
        "SELECT 1 FROM papers WHERE paper_key=?",
        (key,),
    ).fetchone()

    if present and not rescreen_existing:
        with db:
            db.execute(
                "INSERT OR IGNORE INTO paper_peptides VALUES (?,?)",
                (key, peptide),
            )

        return 0, 0, "cached_not_rescreened"

    root = None
    xml_status = "no_fulltext_or_disabled"

    if fulltext and re.fullmatch(r"PMC\d+", pmcid):
        try:
            root = get(
                f"{API}/{pmcid}/fullTextXML",
                as_xml=True,
            )
            xml_status = "retrieved"

        except HTTPError as exc:
            xml_status = "unavailable_http_" + str(exc.code)

        except (
            URLError,
            TimeoutError,
            ET.ParseError,
            ValueError,
        ) as exc:
            xml_status = "fetch_error_" + type(exc).__name__

    kind, reason = classify_publication(p, root)

    blocks = (
        content_blocks(root)
        if root is not None and kind == "original_candidate"
        else []
    )

    new_count = 0

    with db:
        db.execute(
            "INSERT OR IGNORE INTO papers "
            "VALUES (?,?,?,?,?,?,?,?,?,?)",
            (
                key,
                clean(p.get("title")),
                clean(p.get("pubYear")),
                doi,
                clean(p.get("pmid")),
                pmcid,
                url,
                abstract,
                xml_status,
                now(),
            ),
        )

        db.execute(
            "INSERT OR IGNORE INTO paper_peptides VALUES (?,?)",
            (key, peptide),
        )

        db.execute(
            """
            INSERT INTO paper_screening_v2
            VALUES (?,?,?,?,?)
            ON CONFLICT(paper_key) DO UPDATE SET
                study_type=excluded.study_type,
                classification_reason=excluded.classification_reason,
                metadata_json=excluded.metadata_json,
                assessment_at=excluded.assessment_at
            """,
            (
                key,
                kind,
                reason,
                json.dumps(p, ensure_ascii=False),
                now(),
            ),
        )

        for section, paragraph, previous in blocks:
            items = make_candidates(
                paragraph,
                previous,
                peptide,
                section,
            )

            for sentence, code, payload, cited in items:
                payload.update(
                    Reference=clean(p.get("title")),
                    DOI=doi or NA,
                    Paper_URL=url,
                    Study_Year=clean(p.get("pubYear")),
                )

                seed = "\x1f".join(
                    [key, peptide, code, "results_v2", sentence]
                )

                candidate_id = hashlib.sha256(
                    seed.encode()
                ).hexdigest()[:24]

                cursor = db.execute(
                    """
                    INSERT OR IGNORE INTO extraction_candidates
                    VALUES (?,?,?,?,?,?,?,?)
                    """,
                    (
                        candidate_id,
                        key,
                        peptide,
                        "results_v2",
                        sentence,
                        json.dumps(
                            payload,
                            ensure_ascii=False,
                        ),
                        "UNVERIFIED_EXTRACTION",
                        now(),
                    ),
                )

                new_count += int(cursor.rowcount > 0)

                db.execute(
                    """
                    INSERT OR IGNORE INTO candidate_context_v2
                    VALUES (?,?,?,?,?,?,?,?,?,?)
                    """,
                    (
                        candidate_id,
                        key,
                        peptide,
                        section,
                        paragraph,
                        previous,
                        "explicit_same_sentence",
                        cited,
                        int(not cited),
                        now(),
                    ),
                )

    return (
        int(not present),
        new_count,
        kind + ":" + xml_status,
    )


def collect(
    db,
    peptide,
    limit,
    fulltext=True,
    rescreen_existing=False,
):
    query = query_for(peptide)

    data = get(
        API + "/search",
        {
            "query": query,
            "format": "json",
            "resultType": "core",
            "pageSize": min(limit, 100),
            "cursorMark": "*",
        },
    )

    found = (
        data.get("resultList", {})
        .get("result", [])
    )[:limit]

    new_papers = 0
    new_candidates = 0
    statuses = {}

    for paper in found:
        n, c, status = save_paper(
            db,
            paper,
            peptide,
            fulltext,
            rescreen_existing,
        )

        new_papers += n
        new_candidates += c
        statuses[status] = statuses.get(status, 0) + 1

        time.sleep(0.15)

    print(
        json.dumps(
            {
                "peptide": peptide,
                "api_total_matches": data.get("hitCount"),
                "examined": len(found),
                "new_papers": new_papers,
                "new_strict_results_candidates": new_candidates,
                "statuses": statuses,
            },
            ensure_ascii=False,
        )
    )


def report(db):
    tables = (
        "excel_screening",
        "excel_recommendations",
        "papers",
        "extraction_candidates",
        "paper_screening_v2",
        "candidate_context_v2",
    )

    for table in tables:
        count = db.execute(
            "SELECT COUNT(*) FROM " + table
        ).fetchone()[0]

        print(table, count)

    types = db.execute(
        """
        SELECT study_type, COUNT(*)
        FROM paper_screening_v2
        GROUP BY study_type
        """
    ).fetchall()

    print("paper_types:", types)
    print(
        "Note: old candidates remain intact; "
        "v2 candidates are not validated efficacy evidence."
    )


def self_test():
    primary = {
        "title": "Experimental hydrogel for wound healing",
        "pubTypeList": {"pubType": ["Journal Article"]},
    }

    review = {
        "title": "Regenerative peptides: a review",
        "pubTypeList": {"pubType": ["Review"]},
    }

    assert classify_publication(review)[0] == "review"

    xml = ET.fromstring("""
    <article article-type="research-article">
      <body>
        <sec>
          <title>Methods</title>
          <p>GHK increased wound closure.</p>
        </sec>
        <sec>
          <title>Results</title>
          <p>GHK treatment increased wound closure compared with control.</p>
          <p>The scaffold improved migration.</p>
        </sec>
      </body>
    </article>
    """)

    assert classify_publication(
        primary, xml
    )[0] == "original_candidate"

    blocks = content_blocks(xml)

    assert len(blocks) == 2

    assert len(
        make_candidates(
            blocks[0][1], "", "GHK", blocks[0][0]
        )
    ) == 1

    assert len(
        make_candidates(
            blocks[1][1], "", "GHK", blocks[1][0]
        )
    ) == 0

    discussion = ET.fromstring("""
    <article>
      <body>
        <sec>
          <title>Discussion</title>
          <p>GHK improved healing.</p>
        </sec>
      </body>
    </article>
    """)

    assert not content_blocks(discussion)

    print(
        "SELF-TEST PASS: review filter, "
        "Results isolation, explicit peptide, "
        "no false attribution"
    )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--peptide",
        choices=sorted(ALIASES),
        default="GHK",
    )
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--limit", type=int, default=10)
    parser.add_argument("--no-fulltext", action="store_true")

    parser.add_argument(
        "--rescreen-existing",
        action="store_true",
        help="Reinspect prior papers; preserve old evidence candidates",
    )

    parser.add_argument("--report", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--db", type=Path, default=DB)
    parser.add_argument("--excel", type=Path, default=EXCEL)

    args = parser.parse_args()

    if args.self_test:
        self_test()
        return

    if not 1 <= args.limit <= 50:
        parser.error("--limit must be 1..50")

    if not args.excel.is_file():
        parser.error("Original Excel missing: " + str(args.excel))

    db = sqlite3.connect(args.db, timeout=30)

    try:
        initialize(db)

        print("IMPORTED", import_excel(db, args.excel))

        if not args.report:
            peptides = (
                ALIASES if args.all else [args.peptide]
            )

            for peptide in peptides:
                collect(
                    db,
                    peptide,
                    args.limit,
                    not args.no_fulltext,
                    args.rescreen_existing,
                )

        report(db)

    finally:
        db.close()


if __name__ == "__main__":
    main()
