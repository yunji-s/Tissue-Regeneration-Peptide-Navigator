"""Tissue Regeneration Peptide Navigator — read-only evidence dashboard.

Deploy main file: streamlit_app_v2.py
Data build remains separate: python build_database.py --peptide GHK --limit 10
The Excel is never edited. Unverified automatic extractions are not ranked.
"""
from __future__ import annotations

import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

import pandas as pd
import streamlit as st

ROOT = Path(__file__).resolve().parent
DATABASE = ROOT / "navigator.sqlite3"
PEPTIDES = {
    "GHK": "P001", "RADA16-I": "P002", "QK": "P003",
    "IKVAV": "P004", "KLD-12": "P005", "RGD": "P006", "YIGSR": "P007",
}
FUNCTIONS = {"WOUND": "Wound healing", "ANGIO": "Angiogenesis", "NERVE": "Nerve regeneration"}
MISSING = {"", "not reported", "확인필요", "확인 필요", "nan", "none", "null", "n/a"}


def value(raw):
    if raw is None:
        return ""
    return str(raw).strip()


def meaningful(raw):
    return value(raw).lower() not in MISSING


def shown(raw):
    return value(raw) if meaningful(raw) else "확인 필요"


def doi_id(raw):
    return re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value(raw), flags=re.I).lower().strip()


def identity(record):
    """Use DOI first; fall back to study/reference ID. Do not count rows as studies."""
    doi = doi_id(record.get("DOI"))
    if doi:
        return "DOI:" + doi
    study = value(record.get("Study_ID") or record.get("Reference_ID"))
    if study:
        return "STUDY:" + study
    return ""  # Unknown provenance is not an independent study.


def paper_url(record):
    url = value(record.get("Paper_URL"))
    if url.startswith(("https://", "http://")):
        return url
    doi = doi_id(record.get("DOI"))
    return "https://doi.org/" + doi if doi else ""


def classification(record):
    status = value(record.get("Evidence_Status")).lower()
    eligibility = value(record.get("Recommendation_Eligibility")).lower()
    direction = value(record.get("Observed_Direction") or record.get("Effect_Direction")).lower()
    if "non-significant" in status or "nonsignificant" in status or "no significant" in direction:
        return "효과 개선 미입증", "통계적으로 유의한 개선이 보고되지 않은 기록"
    if "unverified" in status or "ineligible" in eligibility or "context only" in eligibility:
        return "추천 제외 · 맥락", "원문 또는 직접 효능 근거가 부족하거나 추천 대상이 아닌 기록"
    if "pending" in status or "caveat" in eligibility or "partial" in status:
        return "검증 필요 후보", "실험 조건·비교군·통계 또는 정량값의 재확인 필요"
    return "검토 후보", "기존 DB 기록으로, 독립적인 원문 재검증은 완료되지 않음"


def join_key(record):
    return (value(record.get("Variant_ID")), doi_id(record.get("DOI")))


def attach_screening(rec, screening_rows):
    """Match variant ID + DOI, then variant name + DOI; never mix different papers."""
    exact = [s for s in screening_rows if meaningful(rec.get("Variant_ID"))
             and join_key(rec) == join_key(s)]
    if exact:
        return exact
    return [s for s in screening_rows
            if value(s.get("Modified_Peptide")).casefold() == value(rec.get("Variant_Name")).casefold()
            and doi_id(s.get("DOI")) == doi_id(rec.get("DOI"))]


def variant_key(rec):
    name = value(rec.get("Variant_Name"))
    return (name.casefold(), value(rec.get("Modification_Method")).casefold(),
            value(rec.get("Formulation_Method")).casefold())


def display_sequence(rec, screenings):
    if screenings:
        sequences = [value(x.get("Modified_Sequence")) for x in screenings if meaningful(x.get("Modified_Sequence"))]
        if sequences:
            seq = sequences[0]
            if any(word in seq.lower() for word in ("검증 필요", "확인필요", "별도 표시")):
                return seq + "  ·  정확한 화학적 구조/서열 검증 필요"
            return seq
    return "변형체의 정확한 서열/구조: 확인 필요"


def reliable_study_count(recs):
    return len({identity(r) for r in recs if identity(r)})


def summary_assessment(recs):
    categories = [classification(r)[0] for r in recs]
    if all(c in {"효과 개선 미입증", "추천 제외 · 맥락"} for c in categories):
        return "근거 부족"
    n = reliable_study_count(recs)
    if n < 2:
        return "단일 연구 근거 · 독립 재현 필요" if n == 1 else "연구 식별정보 부족"
    # Multi-paper evidence does not automatically mean concordant or verified evidence.
    return f"독립 문헌 {n}편 · 결과 일관성 검증 필요"


@st.cache_data(show_spinner=False)
def load_database(db_path, stamp, peptide, function_code):
    del stamp
    with sqlite3.connect(db_path) as conn:
        recs = [json.loads(r[0]) for r in conn.execute(
            "SELECT payload_json FROM excel_recommendations WHERE peptide_id=? AND function_code=?",
            (PEPTIDES[peptide], function_code))]
        screening = [json.loads(r[0]) for r in conn.execute(
            "SELECT payload_json FROM excel_screening WHERE peptide=? AND function_code=?",
            (peptide, function_code))]
        rows = conn.execute(
            "SELECT screening_json FROM extraction_candidates WHERE peptide=?", (peptide,)).fetchall()
        automatic = []
        for (payload,) in rows:
            record = json.loads(payload)
            if value(record.get("Target_Function_Code")) == function_code:
                automatic.append(record)
        total_papers = conn.execute("SELECT COUNT(*) FROM papers").fetchone()[0]
    return recs, screening, automatic, total_papers


def show_paper(record, index=None):
    title = value(record.get("Reference")) or "원저 제목 확인 필요"
    year = value(record.get("Year") or record.get("Study_Year"))
    prefix = f"{index}. " if index is not None else ""
    st.markdown(f"**{prefix}{title}**" + (f" ({year})" if year else ""))
    doi = doi_id(record.get("DOI"))
    url = paper_url(record)
    if doi:
        st.caption(f"DOI: {doi}")
    if url:
        st.link_button("원문 / DOI 열기", url, key="paper_" + str(index) + "_" + re.sub(r"\W", "", doi or url)[-45:])
    else:
        st.caption("원문 링크: 확인 필요")


st.set_page_config(page_title="Tissue Regeneration Peptide Navigator", page_icon="🧬", layout="wide")
st.title("🧬 Tissue Regeneration Peptide Navigator")
st.caption("기능별 후보를 기존 3차 DB에서 조회합니다. 자동 추출 문장은 별도 검증 전까지 추천에 반영하지 않습니다.")

if not DATABASE.is_file():
    st.error("navigator.sqlite3가 없습니다. Codespaces에서 DB 생성 후 GitHub에 commit/push했는지 확인하십시오.")
    st.stop()

c1, c2 = st.columns(2)
with c1:
    peptide = st.selectbox("기준 펩타이드", list(PEPTIDES))
with c2:
    function = st.selectbox("목표 기능", list(FUNCTIONS), format_func=lambda x: FUNCTIONS[x])

if st.button("기능별 후보 조회", use_container_width=True, type="primary"):
    try:
        recs, screening, automatic, paper_total = load_database(
            str(DATABASE), DATABASE.stat().st_mtime_ns, peptide, function)
    except (sqlite3.Error, ValueError, KeyError, json.JSONDecodeError) as exc:
        st.error(f"DB 읽기 실패: {exc}")
        st.stop()

    st.caption(f"기존 추천 {len(recs)}행 · 기존 Screening {len(screening)}행 · 자동 추출 {len(automatic)}문장(미검증) · 전체 수집 논문 {paper_total}편")
    st.subheader("추천 검토 후보")
    st.info("표시 순서는 효능 순위가 아닙니다. 검증 완료된 독립 연구가 부족한 후보는 '우선 검토 후보'로만 표시합니다.")

    groups = defaultdict(list)
    for record in recs:
        groups[variant_key(record)].append(record)

    candidate_groups = []
    excluded_groups = []
    for key, group in groups.items():
        usable = [r for r in group if classification(r)[0] in {"검증 필요 후보", "검토 후보"}]
        if usable:
            candidate_groups.append((key, usable, group))
        else:
            excluded_groups.append((key, group))

    # Stable presentation only; never pretend this order is effect-size ranking.
    candidate_groups.sort(key=lambda entry: (-reliable_study_count(entry[1]), entry[0][0]))
    if not candidate_groups:
        st.warning("현재 목표 기능에 표시할 검토 후보가 없습니다. 근거 없이 3개를 채우지 않습니다.")

    for index, (key, usable, all_records) in enumerate(candidate_groups[:3], start=1):
        representative = usable[0]
        connected_screening = []
        for record in all_records:
            connected_screening.extend(attach_screening(record, screening))
        unique_screening = {value(r.get("Record_ID")): r for r in connected_screening}
        connected_screening = list(unique_screening.values())
        name = shown(representative.get("Variant_Name"))
        with st.container(border=True):
            st.markdown(f"### 검토 후보 {index} · {name}")
            st.caption(summary_assessment(usable))
            a, b = st.columns([1.3, 1])
            with a:
                st.markdown("**변형 서열/구조**")
                st.write(display_sequence(representative, connected_screening))
                st.markdown("**Modification / Formulation**")
                st.write(shown(representative.get("Modification_Method")))
                st.write(shown(representative.get("Formulation_Method")))
            with b:
                st.markdown("**문헌에서 보고된 기능**")
                st.write(shown(representative.get("Target_Function_EN")) or FUNCTIONS[function])
                st.markdown("**평가 지표**")
                st.write(" · ".join(sorted({shown(r.get("Outcome_Metric")) for r in usable})))
                st.markdown("**독립 원저 수**")
                st.write(f"{reliable_study_count(usable)}편 (DOI / Study ID 기준, 동일 논문 중복 제외)")
            st.markdown("**대표 원문**")
            url = paper_url(representative)
            st.write(value(representative.get("Reference")) or "서지정보 확인 필요")
            if url:
                st.link_button("대표 원문 확인", url, key=f"top_{index}")
            else:
                st.caption("대표 원문 링크 확인 필요")
            with st.expander("더보기 · 실험 근거 / 비교군 / 여러 원문 / 한계"):
                seen_papers = set()
                for n, record in enumerate(usable, start=1):
                    st.markdown(f"**근거 기록 {n}**")
                    st.write("효과 방향:", shown(record.get("Observed_Direction")))
                    st.write("비교군:", shown(record.get("Comparator")))
                    st.write("평가 방법/지표:", shown(record.get("Outcome_Metric")))
                    st.write("근거 요약:", shown(record.get("Evidence_Detail") or record.get("Short_Summary")))
                    st.write("통계/검증 상태:", shown(record.get("Evidence_Status")))
                    st.write("한계:", shown(record.get("Limitations")))
                    ident = identity(record)
                    if ident and ident not in seen_papers:
                        seen_papers.add(ident)
                        st.markdown("원문: " + (f"[링크]({paper_url(record)})" if paper_url(record) else "확인 필요"))
                    st.divider()
                if connected_screening:
                    st.markdown("**기존 Screening DB 원자료**")
                    cols = ["Modified_Peptide", "Modified_Sequence", "Evaluation_Method", "Quantitative_Result",
                            "Statistical_Significance", "Comparator", "DOI", "Source_Verification"]
                    st.dataframe(pd.DataFrame(connected_screening).reindex(columns=cols), use_container_width=True, hide_index=True)
                if len(seen_papers) < 2:
                    st.warning("독립 원저가 2편 미만입니다. 단일 논문의 여러 실험 기록을 반복 재현 근거로 계산하지 않습니다.")

    with st.expander("추천에 사용하지 않은 근거 · 제외 사유와 다음 검증 작업", expanded=False):
        st.write("비유의 결과, 직접 효능이 불명확한 결과, 비교군이 섞인 결과는 삭제하지 않고 분리 보관합니다.")
        excluded_rows = []
        for _, group in excluded_groups:
            for record in group:
                c, reason = classification(record)
                excluded_rows.append({
                    "변형체": shown(record.get("Variant_Name")), "분류": c, "사유": reason,
                    "관찰 결과": shown(record.get("Observed_Direction")),
                    "DOI": doi_id(record.get("DOI")),
                    "다음 작업": "원저·비교군·통계값 검토 후 상태 갱신",
                })
        if excluded_rows:
            st.dataframe(pd.DataFrame(excluded_rows), use_container_width=True, hide_index=True)
        else:
            st.caption("현재 조건에서 제외된 기존 추천 기록은 없습니다.")
        st.markdown("**자동 추출된 미검증 근거**")
        st.write(f"{len(automatic)}문장 · 모두 추천 및 독립 원저 수 계산에서 제외")
        if automatic:
            columns = ["Original_Peptide", "Candidate_Statement", "Source_Verification", "DOI", "Evidence_Status"]
            table = pd.DataFrame(automatic).reindex(columns=columns)
            st.dataframe(table, use_container_width=True, hide_index=True)
            st.download_button("미검증 근거 CSV 다운로드", table.to_csv(index=False).encode("utf-8-sig"),
                               "unverified_evidence.csv", "text/csv")
        st.caption("검증 작업: (1) 해당 peptide intervention 확인 → (2) 결과 지표 및 비교군 연결 → (3) 정량값/p-value 확인 → (4) DOI 중복 검사 → (5) 유효 근거로 승격")

    if len(candidate_groups) > 3:
        st.caption(f"추가 검토 후보 {len(candidate_groups)-3}개는 이번 화면에서 표시하지 않았습니다. 화면의 3개는 임시 효능 순위가 아닙니다.")
