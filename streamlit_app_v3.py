"""Tissue Regeneration Peptide Navigator — literature modification library.

Read-only viewer of the EXISTING 3rd screening database.
No ranking, no new paper retrieval, no automatic claim verification.
Joins v13 evidence by Record_ID, but never treats an automated label as human-verified.
Run: python -m streamlit run streamlit_app_v3.py
"""
from __future__ import annotations

import io
import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
DB = ROOT / "navigator.sqlite3"
PEPTIDES = ("GHK", "RADA16-I", "QK", "IKVAV", "KLD-12", "RGD", "YIGSR")
EMPTY = {"", "none", "nan", "not reported", "확인필요", "확인 필요", "-"}
FUNCTION_NAMES = {
    "WOUND": "상처 치유", "ANGIO": "혈관신생", "NERVE": "신경재생",
    "CARTILAGE": "연골재생", "ADHESION": "세포 부착", "MIGRATE": "세포 이동",
    "PROLIF": "세포 증식", "ECMREMOD": "ECM 재형성", "STABILITY": "안정성",
    "MYOGEN": "근육 형성", "APOPTOSIS": "세포사멸", "DIFF": "세포 분화",
    "DELIVERY": "전달", "INFLAM": "염증", "ANTIINFLAM": "항염증",
}


def txt(value: object) -> str:
    return str(value if value is not None else "").strip()


def visible(value: object) -> str:
    text = txt(value)
    return "미확인" if text.lower() in EMPTY else text


def doi_normalized(value: object) -> str:
    return re.sub(r"^https?://(?:dx\.)?doi\.org/", "", txt(value), flags=re.I).lower()


def reference_url(row: dict) -> str:
    url = txt(row.get("Paper_URL"))
    if url.startswith(("https://", "http://")):
        return url
    doi = doi_normalized(row.get("DOI"))
    return f"https://doi.org/{doi}" if doi else ""


def included(parent: str, selected: str) -> bool:
    if selected == "7종 전체":
        return True
    # Composite RGD + YIGSR retains its original name and is clearly marked.
    return selected in [part.strip() for part in re.split(r"\s*\+\s*", parent)]


@st.cache_data(show_spinner=False)
def read_screening(path: str, mtime_ns: int, size: int) -> list[dict]:
    uri = Path(path).resolve().as_uri() + "?mode=ro"
    with sqlite3.connect(uri, uri=True) as conn:
        found = conn.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='excel_screening'"
        ).fetchone()
        if not found:
            raise RuntimeError("기존 DB에 excel_screening 테이블이 없습니다.")
        items = []
        for (raw,) in conn.execute("SELECT payload_json FROM excel_screening"):
            try:
                item = json.loads(raw)
                if isinstance(item, dict):
                    items.append(item)
            except (ValueError, TypeError):
                continue
    return items


@st.cache_data(show_spinner=False)
def read_evidence(path: str, mtime_ns: int, size: int) -> tuple[dict, dict]:
    """Read separate v13 outputs, if available; no database writes."""
    uri = Path(path).resolve().as_uri() + "?mode=ro"
    summaries, candidates = {}, {}
    with sqlite3.connect(uri, uri=True) as conn:
        conn.row_factory = sqlite3.Row
        tables = {row[0] for row in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "qualitative_summaries_v2" in tables:
            for item in conn.execute("SELECT * FROM qualitative_summaries_v2"):
                data = dict(item)
                rid = txt(data.get("record_id"))
                if rid:
                    summaries[rid] = data
        if "qualitative_evidence_candidates_v1" in tables:
            for item in conn.execute("SELECT * FROM qualitative_evidence_candidates_v1 ORDER BY record_id, rank"):
                data = dict(item)
                rid = txt(data.get("record_id"))
                if rid:
                    candidates.setdefault(rid, []).append(data)
    return summaries, candidates


def evidence_level(record: dict, entry: dict | None) -> str:
    """Presentation-only labels. Do not conflate automated checks with verification."""
    if not entry:
        return "근거 미연결"
    original_doi = doi_normalized(record.get("DOI"))
    evidence_doi = doi_normalized(entry.get("doi"))
    if original_doi and evidence_doi and original_doi != evidence_doi:
        return "검토 필요"
    status = txt(entry.get("verification_status")).upper()
    check = txt(entry.get("certainty_check")).upper()
    quote = txt(entry.get("evidence_quote"))
    summary = txt(entry.get("summary_ko"))
    if not quote or not summary or "FALLBACK" in status or "FLAGGED" in check or "NO_EVIDENCE" in check:
        return "검토 필요"
    if "MODIFICATION_SPECIFIC" in status:
        return "변형체 직접 언급"
    if "STUDY_LEVEL_ONLY" in status:
        return "논문 수준 근거"
    return "검토 필요"


def render_automated_evidence(record: dict, summaries: dict, candidates: dict) -> None:
    rid = txt(record.get("Record_ID"))
    entry = summaries.get(rid)
    level = evidence_level(record, entry)
    st.markdown("**v13 자동 분석 결과 (사람의 원문 검증 전)**")
    if not entry:
        st.caption("분석 결과 없음: 기존 Screening DB 항목은 그대로 보존됩니다.")
    else:
        st.caption(f"근거 구분: {level} · 분석 상태: {visible(entry.get('verification_status'))} · "
                   f"근거 유형: {visible(entry.get('evidence_type'))} · 자동 검사: {visible(entry.get('certainty_check'))}")
        summary = txt(entry.get("summary_ko"))
        quote = txt(entry.get("evidence_quote"))
        if level == "검토 필요":
            st.warning("자동 요약 또는 원문 연결 검토가 필요합니다. 아래 내용을 확정된 효능으로 해석하지 마십시오.")
        if summary:
            st.write("**Research outcome (English, automated):**", summary)
        if quote:
            st.write("**Evidence quote (stored abstract extract):**", quote)
        if txt(entry.get("review_note")):
            st.caption("검토 메모: " + txt(entry["review_note"]))
        if txt(entry.get("source_kind")):
            st.caption("출처 유형: " + txt(entry["source_kind"]) + " · 모델: " + visible(entry.get("model")))
    matches = candidates.get(rid, [])
    if matches:
        with st.expander(f"근거 후보 {len(matches)}개 확인", expanded=False):
            for item in matches:
                st.caption(f"우선순위 {visible(item.get('rank'))} · 선택 {bool(item.get('selected'))} · "
                           f"범위 {visible(item.get('scope'))} · 유형 {visible(item.get('evidence_type'))}")
                st.write(txt(item.get("evidence_quote")))


def export_connected_csv(rows: list[dict], summaries: dict) -> bytes:
    out = []
    for record in rows:
        item = dict(record)
        entry = summaries.get(txt(record.get("Record_ID")))
        item["AI_Evidence_Level"] = evidence_level(record, entry)
        for field, label in (
            ("summary_ko", "AI_Research_Outcome_English"),
            ("evidence_quote", "AI_Evidence_Quote"),
            ("verification_status", "AI_Verification_Status"),
            ("evidence_type", "AI_Evidence_Type"),
            ("certainty_check", "AI_Certainty_Check"),
            ("review_note", "AI_Review_Note"),
        ):
            item[label] = txt(entry.get(field)) if entry else ""
        out.append(item)
    return export_csv(out)


def group_library(rows: list[dict]) -> list[dict]:
    groups: dict[tuple[str, ...], list[dict]] = defaultdict(list)
    for row in rows:
        key = tuple(txt(row.get(column)) for column in (
            "Original_Peptide", "Modified_Peptide", "Modification_Method",
            "Formulation", "Intervention_Class"
        ))
        groups[key].append(row)

    output = []
    for (parent, variant, modification, formulation, category), items in groups.items():
        codes = sorted({txt(r.get("Target_Function_Code")) for r in items if txt(r.get("Target_Function_Code"))})
        dois = sorted({doi_normalized(r.get("DOI")) for r in items if doi_normalized(r.get("DOI"))})
        output.append({
            "parent": parent, "variant": variant, "modification": modification,
            "formulation": formulation, "category": category,
            "functions": codes, "dois": dois, "rows": items,
            "is_composite": "+" in parent,
        })
    # Alphabetical display, deliberately NOT efficacy ranking.
    return sorted(output, key=lambda x: (
        x["parent"].casefold(), x["variant"].casefold(), x["formulation"].casefold()
    ))


def export_csv(rows: list[dict]) -> bytes:
    if not rows:
        return b""
    columns = list(dict.fromkeys(k for row in rows for k in row))
    df = pd.DataFrame([{k: row.get(k, "") for k in columns} for row in rows])
    return df.to_csv(index=False).encode("utf-8-sig")


def render_record(row: dict, summaries: dict, candidates: dict) -> None:
    title = txt(row.get("Reference")) or "출처 미기재"
    year = txt(row.get("Study_Year"))
    st.markdown(f"**{title}{f' ({year})' if year else ''}**")
    st.caption(
        f"Record ID: {visible(row.get('Record_ID'))} · "
        f"기존 DB 근거 상태: {visible(row.get('Evidence_Status'))} · "
        f"원문 검증 메모: {visible(row.get('Source_Verification'))}"
    )
    st.write("**목표 기능:**", visible(row.get("Target_Function")))
    st.write("**평가 방법:**", visible(row.get("Evaluation_Method")))
    st.write("**평가 지표:**", visible(row.get("Measured_Outcome")))
    st.write("**효과 방향:**", visible(row.get("Effect_Direction")))
    st.write("**정량·정성 결과:**", visible(row.get("Quantitative_Result")))
    st.write("**비교군:**", visible(row.get("Comparator")))
    st.write("**통계:**", visible(row.get("Statistical_Significance")))
    if txt(row.get("Effect_Interpretation")).lower() not in EMPTY:
        st.write("**기존 DB 해석:**", row["Effect_Interpretation"])
    if txt(row.get("Reviewer_Notes")).lower() not in EMPTY:
        st.info("검토 메모: " + txt(row["Reviewer_Notes"]))
    render_automated_evidence(row, summaries, candidates)
    url = reference_url(row)
    if url:
        st.link_button("논문 출처 열기", url)
    else:
        st.warning("DOI 또는 출처 링크가 없습니다.")


st.set_page_config(page_title="Peptide Modification Library", page_icon="🧬", layout="wide")
st.title("🧬 Peptide Modification Library")
st.caption("기존 3차 DB 기반 · 7종 펩타이드 · 변형체 및 제형 문헌 라이브러리")
st.info(
    "현재 화면은 **추천 순위가 아닌 문헌 기록 목록**입니다. "
    "기존 DB와 v13 자동 분석 결과를 구분하여 보여주며, "
    "새 논문 검색·AI 검증·정량 재계산을 수행하지 않습니다."
)

if not DB.is_file():
    st.error("navigator.sqlite3 파일이 없습니다. 기존 프로젝트의 DB 파일을 이 코드와 같은 폴더에 넣으십시오.")
    st.stop()

try:
    info = DB.stat()
    source_rows = read_screening(str(DB), info.st_mtime_ns, info.st_size)
    summaries, candidates = read_evidence(str(DB), info.st_mtime_ns, info.st_size)
except Exception as exc:
    st.error(f"기존 DB 조회 실패: {exc}")
    st.stop()

if not source_rows:
    st.warning("Screening DB에 조회 가능한 기록이 없습니다.")
    st.stop()

library = group_library(source_rows)
all_functions = sorted({code for group in library for code in group["functions"]})
all_categories = sorted({group["category"] for group in library if group["category"]})

with st.sidebar:
    st.header("라이브러리 필터")
    parent = st.selectbox("원형 펩타이드", ["7종 전체", *PEPTIDES])
    selected_function = st.selectbox(
        "목표 기능", ["전체 기능", *all_functions],
        format_func=lambda x: x if x == "전체 기능" else f"{FUNCTION_NAMES.get(x, x)} ({x})"
    )
    selected_category = st.selectbox("연구 유형", ["전체 유형", *all_categories])
    keyword = st.text_input("변형체·소재·지표 검색", placeholder="예: GHK-Cu, hydrogel, CD31")
    verified_only = st.checkbox("기존 DB에 'Verified'로 표기된 기록만", value=False)
    st.caption("이 체크박스는 기존 DB의 표시를 기준으로 하며 원저 재검증을 의미하지 않습니다.")
    selected_level = st.selectbox("AI 근거 수준", ["전체", "변형체 직접 언급", "논문 수준 근거", "검토 필요", "근거 미연결"])

filtered = []
for group in library:
    if not included(group["parent"], parent):
        continue
    if selected_function != "전체 기능" and selected_function not in group["functions"]:
        continue
    if selected_category != "전체 유형" and selected_category != group["category"]:
        continue
    records = group["rows"]
    if selected_function != "전체 기능":
        records = [r for r in records if txt(r.get("Target_Function_Code")) == selected_function]
    if verified_only:
        records = [r for r in records if "verified" in txt(r.get("Evidence_Status")).lower()
                   and "unverified" not in txt(r.get("Evidence_Status")).lower()]
    if selected_level != "전체":
        records = [r for r in records if evidence_level(r, summaries.get(txt(r.get("Record_ID")))) == selected_level]
    if not records:
        continue
    if keyword:
        blob = " ".join(txt(v) for row in records for v in row.values())
        blob += " " + " ".join([group["variant"], group["modification"], group["formulation"]])
        if keyword.casefold() not in blob.casefold():
            continue
    filtered.append({**group, "rows": records})

matched_records = [r for group in filtered for r in group["rows"]]
unique_dois = {doi_normalized(row.get("DOI")) for row in matched_records if doi_normalized(row.get("DOI"))}

c1, c2, c3 = st.columns(3)
c1.metric("표시된 변형체·제형 조합", len(filtered))
c2.metric("연결된 Screening 기록", len(matched_records))
c3.metric("고유 DOI", len(unique_dois))
st.caption("논문 한 편에 여러 결과가 있을 수 있어 DOI 수와 실험 기록 수는 서로 다릅니다.")

st.download_button(
    "📥 현재 검색 결과 CSV 다운로드",
    data=export_connected_csv(matched_records, summaries),
    file_name="peptide_modification_library_filtered.csv",
    mime="text/csv", disabled=not matched_records
)

st.subheader("문헌에 보고된 Modification·Formulation")
st.caption("원형 펩타이드·변형체·제형명에 따른 가나다/알파벳 순서이며, 효능의 우열을 나타내지 않습니다.")

if not filtered:
    st.info("조건에 해당하는 기존 DB 기록이 없습니다. 필터를 조정해 보십시오.")

for idx, group in enumerate(filtered):
    with st.container(border=True):
        name = visible(group["variant"])
        st.markdown(f"#### {name}")
        sub = " + ".join([group["parent"], group["category"]])
        st.caption(sub + (" · 복합 펩타이드 기록" if group["is_composite"] else ""))
        left, right = st.columns(2)
        with left:
            st.write("**Modification 방법:**", visible(group["modification"]))
            st.write("**Formulation:**", visible(group["formulation"]))
        with right:
            orig = next((r.get("Original_Sequence") for r in group["rows"] if txt(r.get("Original_Sequence"))), "")
            changed = next((r.get("Modified_Sequence") for r in group["rows"] if txt(r.get("Modified_Sequence"))), "")
            st.write("**원형 서열:**", visible(orig))
            st.write("**변형 서열/구조:**", visible(changed))
        st.write("**연구된 기능:**", ", ".join(FUNCTION_NAMES.get(code, code) for code in group["functions"]))
        st.caption(f"연결된 기록 {len(group['rows'])}건 · 고유 DOI {len({doi_normalized(r.get('DOI')) for r in group['rows'] if doi_normalized(r.get('DOI'))})}개")
        with st.expander("논문별 정량·정성 근거 및 원문 출처"):
            for i, record in enumerate(group["rows"]):
                if i:
                    st.divider()
                render_record(record, summaries, candidates)

st.divider()
st.caption(
    f"현재 DB의 전체 Screening 기록 {len(source_rows)}건 · "
    f"자동 분석 결과 연결 {sum(1 for r in source_rows if txt(r.get('Record_ID')) in summaries)}건 · "
    "AI 자동 분류와 요약은 원저 대조 및 연구자 검토 전의 잠정 자료입니다."
)
