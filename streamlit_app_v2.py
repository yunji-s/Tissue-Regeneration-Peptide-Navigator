"""Tissue Regeneration Peptide Navigator — function-first 7-peptide evidence browser.

SQLite source: build_database.py --all --limit 10
Unverified Europe PMC sentences are never promoted to verified recommendations.
"""
from __future__ import annotations

import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

import pandas as pd
import streamlit as st

DB_PATH = Path(__file__).resolve().parent / "navigator.sqlite3"
PEPTIDES = ("GHK", "RADA16-I", "QK", "IKVAV", "KLD-12", "RGD", "YIGSR")
PEPTIDE_IDS = {name: f"P{i:03d}" for i, name in enumerate(PEPTIDES, 1)}
FUNCTIONS = {"WOUND": "Wound healing", "ANGIO": "Angiogenesis", "NERVE": "Nerve regeneration",
             "CARTILAGE": "Cartilage regeneration", "MIGRATE": "Cell migration", "ADHESION": "Cell adhesion",
             "PROLIF": "Cell proliferation", "INFLAM": "Inflammatory response", "ECMREMOD": "ECM remodeling"}
MISSING = {"", "not reported", "nan", "none", "null", "n/a", "확인필요", "확인 필요", "-"}


def value(x):
    return "" if x is None else str(x).strip()


def valid(x):
    return value(x).lower() not in MISSING


def show(x):
    return value(x) if valid(x) else "확인 필요"


def doi(x):
    return re.sub(r"^https?://(?:dx\.)?doi\.org/", "", value(x), flags=re.I).lower()


def link(row):
    url = value(row.get("Paper_URL") or row.get("paper_url"))
    if url.startswith(("https://", "http://")):
        return url
    return "https://doi.org/" + doi(row.get("DOI") or row.get("doi")) if doi(row.get("DOI") or row.get("doi")) else ""


def paper_id(row):
    d = doi(row.get("DOI") or row.get("doi"))
    if d:
        return "DOI:" + d
    for field in ("PMID", "pmid", "PMCID", "pmcid"):
        if valid(row.get(field)):
            return field.upper() + ":" + value(row[field])
    # No paper ID -> cannot count as a distinct independent publication
    return ""


def kind(row):
    status = value(row.get("Evidence_Status")).lower()
    eligibility = value(row.get("Recommendation_Eligibility")).lower()
    direction = value(row.get("Observed_Direction")).lower()
    if ("non-significant" in status or "nonsignificant" in status
            or "no significant" in direction or "show as non-significant" in eligibility):
        return "negative", "비유의 결과 — 개선 근거로 사용하지 않음"
    if ("unverified" in status or "recheck" in status or "context only" in eligibility
            or "exclude" in eligibility or "ineligible" in eligibility):
        return "excluded", "원문 확인 또는 직접 기능 근거 부족"
    if "verified significant" in status and "eligible (verified)" in eligibility:
        return "verified_db", "기존 DB에서 유의성 확인으로 표시 — 독립 재검증 필요"
    return "limited", "기존 DB의 제한적 근거 — 비교군·정량값 추가 검증 필요"


def variant_key(row):
    # Keep parent and specific modified/formulation identity separate.
    return (value(row.get("Parent_Peptide")).casefold(),
            value(row.get("Variant_Name")).casefold(),
            value(row.get("Modification_Method")).casefold(),
            value(row.get("Formulation_Method")).casefold())


def related_screening(record, screening):
    rid = value(record.get("Variant_ID"))
    d = doi(record.get("DOI"))
    matched = [s for s in screening if value(s.get("Variant_ID")) == rid and rid
               and doi(s.get("DOI")) == d and d]
    if not matched:
        matched = [s for s in screening if value(s.get("Modified_Peptide")).casefold()
                   == value(record.get("Variant_Name")).casefold() and doi(s.get("DOI")) == d and d]
    return matched


def sequence_text(record, screening):
    for s in screening:
        if valid(s.get("Modified_Sequence")):
            return value(s["Modified_Sequence"])
    return "정확한 변형 서열/구조: 확인 필요"


@st.cache_data(show_spinner=False)
def load_all(path, modified_ns, byte_size):
    # Cache refreshed when the SQLite file changes.
    del modified_ns, byte_size
    with sqlite3.connect(path) as con:
        tables = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        required = {"excel_recommendations", "excel_screening", "papers", "extraction_candidates"}
        if not required.issubset(tables):
            raise ValueError("SQLite schema incomplete; run build_database.py --report")
        rec = [json.loads(r[0]) for r in con.execute("SELECT payload_json FROM excel_recommendations")]
        screen = [json.loads(r[0]) for r in con.execute("SELECT payload_json FROM excel_screening")]
        papers = [dict(zip(("paper_key", "title", "year", "doi", "pmid", "pmcid", "paper_url", "fulltext_status"), r))
                  for r in con.execute("SELECT paper_key,title,year,doi,pmid,pmcid,paper_url,fulltext_status FROM papers")]
        links = con.execute("SELECT paper_key,peptide FROM paper_peptides").fetchall()
        candidates = [dict(zip(("paper_key", "peptide", "source_section", "source_text", "screening_json", "evidence_status"), r))
                      for r in con.execute("SELECT paper_key,peptide,source_section,source_text,screening_json,evidence_status FROM extraction_candidates")]
    associated = defaultdict(set)
    for key, peptide in links:
        associated[key].add(peptide)
    for paper in papers:
        paper["peptides"] = associated[paper["paper_key"]]
    return rec, screen, papers, candidates


def category_for_candidate(row):
    try:
        return value(json.loads(row["screening_json"]).get("Target_Function_Code"))
    except (ValueError, TypeError):
        return ""


def display_reference(row):
    name = show(row.get("Reference") or row.get("title"))
    d = doi(row.get("DOI") or row.get("doi"))
    url = link(row)
    if url:
        st.markdown(f"[{name}]({url})" + (f" · DOI `{d}`" if d else ""))
    else:
        st.write(name)


st.set_page_config(page_title="Tissue Regeneration Peptide Navigator", page_icon="🧬", layout="wide")
st.title("🧬 Tissue Regeneration Peptide Navigator")
st.caption("7종 펩타이드 · 기능 기반 Modification / Formulation 검토 후보 · Europe PMC 문헌 통합 조회")

if not DB_PATH.is_file():
    st.error("navigator.sqlite3 파일이 없습니다. Codespaces에서 `python build_database.py --all --limit 10`을 실행한 뒤 GitHub에 DB를 업로드하십시오.")
    st.stop()

try:
    stat = DB_PATH.stat()
    records, screenings, papers, auto = load_all(str(DB_PATH), stat.st_mtime_ns, stat.st_size)
except (sqlite3.Error, ValueError, KeyError, json.JSONDecodeError) as err:
    st.error(f"데이터베이스 읽기 실패: {err}")
    st.stop()

left, right = st.columns([1.5, 1])
with left:
    function_code = st.selectbox("① 원하는 조직재생 기능", list(FUNCTIONS),
                                 format_func=lambda x: FUNCTIONS[x])
with right:
    scope = st.selectbox("② 검색 범위", ["7종 전체", *PEPTIDES])
selected = set(PEPTIDES) if scope == "7종 전체" else {scope}

if st.button("문헌 DB 통합 조회 및 후보 추천", type="primary", use_container_width=True):
    target_ids = {PEPTIDE_IDS[name] for name in selected}
    matched = [r for r in records if value(r.get("Target_Function_Code")) == function_code
               and (value(r.get("Peptide_ID")) in target_ids
                    or value(r.get("Parent_Peptide")) in selected
                    or any(pid in value(r.get("Peptide_ID")).split(";") for pid in target_ids))]
    reference_screening = [s for s in screenings if value(s.get("Target_Function_Code")) == function_code]
    relevant_papers = [p for p in papers if p["peptides"].intersection(selected)]
    relevant_auto = [x for x in auto if x["peptide"] in selected and category_for_candidate(x) == function_code]
    st.caption(f"기존 추천 기록 {len(matched)}건 · 수집된 Europe PMC 원저 후보 {len(relevant_papers)}편(기능 검증 전 포함) · 기능 연관 자동 추출 문장 {len(relevant_auto)}건")

    grouped = defaultdict(list)
    for r in matched:
        grouped[variant_key(r)].append(r)
    ranked, excluded = [], []
    for key, rows in grouped.items():
        positive = [r for r in rows if kind(r)[0] in {"verified_db", "limited"}]
        for r in rows:
            if kind(r)[0] not in {"verified_db", "limited"}:
                excluded.append(r)
        if not positive:
            continue
        identifiers = {paper_id(r) for r in positive if paper_id(r)}
        # Priority is for review, not a quantified biological efficacy score.
        best_label = 1 if any(kind(r)[0] == "verified_db" for r in positive) else 0
        ranked.append((best_label, len(identifiers), key, positive, rows, identifiers))
    ranked.sort(key=lambda x: (-x[0], -x[1], x[2]))

    st.subheader("추천 검토 후보 TOP 1–3")
    st.info("기존 DB에 등록된 변형·소재화 후보를 비교합니다. 표시 순서는 검토 우선순위이며 효능이 가장 높다는 뜻이 아닙니다. Europe PMC 미검증 문장은 순위에 가산하지 않습니다.")
    if not ranked:
        st.warning("현재 조건에서 직접 검토할 수 있는 추천 기록이 부족합니다. 무리하게 후보 3개를 채우지 않습니다.")

    for rank, (_, distinct, _, supporting, all_rows, ids) in enumerate(ranked[:3], 1):
        r = supporting[0]
        parent = show(r.get("Parent_Peptide"))
        name = show(r.get("Variant_Name"))
        specific_screen = {}
        for item in all_rows:
            for s in related_screening(item, reference_screening):
                specific_screen[value(s.get("Record_ID"))] = s
        screen_rows = list(specific_screen.values())
        # Match collected Europe PMC metadata ONLY by exact DOI of this variant's existing records.
        linked = [p for p in relevant_papers if doi(p.get("doi"))
                  and doi(p.get("doi")) in {doi(e.get("DOI")) for e in all_rows}]
        with st.container(border=True):
            st.markdown(f"### {rank}. {name}  ·  원형 {parent}")
            st.caption(f"독립 원저 식별 {distinct}편 · " +
                       ("독립 재현과 결과 일관성 추가 검증 필요" if distinct >= 2 else "단일 원저 또는 문헌 식별정보 부족"))
            a, b = st.columns(2)
            with a:
                st.markdown("**변형 서열 / 구조 표기**")
                st.code(sequence_text(r, screen_rows), language=None)
                st.write("**Modification:**", show(r.get("Modification_Method")))
                st.write("**Formulation:**", show(r.get("Formulation_Method")))
            with b:
                st.markdown("**문헌 기반 기대 기능**")
                st.markdown("#### " + FUNCTIONS[function_code])
                st.write("**대표 평가 지표:**", ", ".join(sorted({show(x.get("Outcome_Metric")) for x in supporting})))
                st.write("**근거 분류:**", "기존 DB 확인 항목 포함" if any(kind(x)[0] == "verified_db" for x in supporting) else "제한적 근거")
            st.markdown("**대표 원문**")
            display_reference(r)
            with st.expander("더보기 — 실험 조건·비교군·정량 결과·전체 원문·한계"):
                for j, x in enumerate(all_rows, 1):
                    st.markdown(f"**기존 DB 실험 기록 {j} · {kind(x)[1]}**")
                    st.write("평가 지표:", show(x.get("Outcome_Metric")))
                    st.write("효과 방향:", show(x.get("Observed_Direction")))
                    st.write("비교군:", show(x.get("Comparator")))
                    st.write("근거 설명:", show(x.get("Evidence_Detail") or x.get("Short_Summary")))
                    st.write("한계:", show(x.get("Limitations")))
                    display_reference(x)
                    st.divider()
                if screen_rows:
                    st.markdown("**기존 Screening DB 원자료 (정량 결과·통계 확인용)**")
                    columns = ["Modified_Peptide", "Modified_Sequence", "Evaluation_Method", "Quantitative_Result",
                               "Statistical_Significance", "Comparator", "Source_Verification", "DOI"]
                    st.dataframe(pd.DataFrame(screen_rows).reindex(columns=columns), hide_index=True, use_container_width=True)
                st.markdown("**Europe PMC 연계 원문 (기존 DB와 DOI가 일치하는 경우)**")
                if linked:
                    for p in linked:
                        display_reference(p)
                        st.caption("Europe PMC 수집 메타데이터 연결 · 독립적인 추가 재현 근거로 중복 계산하지 않음")
                else:
                    st.caption("이번 수집 범위 내에서 DOI가 일치하는 추가 메타데이터 없음")

    st.subheader("Europe PMC 추가 문헌 탐색 결과")
    st.caption("7종 전체(또는 선택한 펩타이드)에서 수집된 논문입니다. 기능에 대한 실제 효능·변형체 일치는 아직 검증되지 않았습니다.")
    search_phrase = st.text_input("수집 논문 제목/DOI 검색", placeholder="예: GHK, hydrogel, angiogenesis", key="paper_query")
    shown_papers = [p for p in relevant_papers if not search_phrase or search_phrase.lower() in (p["title"] + " " + value(p["doi"])).lower()]
    if shown_papers:
        paper_table = [{"펩타이드(검색 연결)": ", ".join(sorted(p["peptides"])), "논문": p["title"],
                        "연도": p["year"], "DOI": p["doi"], "원문 링크": p["paper_url"],
                        "원문 XML": p["fulltext_status"]} for p in shown_papers]
        st.dataframe(pd.DataFrame(paper_table), use_container_width=True, hide_index=True,
                     column_config={"원문 링크": st.column_config.LinkColumn("원문 링크")})
    else:
        st.info("저장된 해당 펩타이드의 신규 논문이 없습니다. Codespaces에서 `python build_database.py --all --limit 10`을 실행하고 SQLite를 GitHub에 push하십시오.")

    with st.expander("추천에 사용하지 않은 근거 / 미검증 추출 문장", expanded=False):
        st.markdown("**기존 DB에서 제외된 결과**")
        if excluded:
            st.dataframe(pd.DataFrame([{"원형": x.get("Parent_Peptide"), "변형체": x.get("Variant_Name"),
                                       "사유": kind(x)[1], "관찰 결과": x.get("Observed_Direction"),
                                       "DOI": x.get("DOI")} for x in excluded]),
                         hide_index=True, use_container_width=True)
        else:
            st.caption("현재 조건에서 분류된 제외 기록 없음")
        st.markdown(f"**Europe PMC 자동 추출 문장 — {len(relevant_auto)}건, 전부 미검증**")
        if relevant_auto:
            rows = []
            for x in relevant_auto:
                info = json.loads(x["screening_json"])
                rows.append({"펩타이드": x["peptide"], "기능": info.get("Target_Function_Code"),
                             "문장": x["source_text"], "출처": x["source_section"],
                             "DOI": info.get("DOI"), "상태": x["evidence_status"]})
            frame = pd.DataFrame(rows)
            st.dataframe(frame, hide_index=True, use_container_width=True)
            st.download_button("미검증 문장 CSV 받기", frame.to_csv(index=False).encode("utf-8-sig"),
                               file_name="unverified_extraction.csv", mime="text/csv")
        else:
            st.caption("선택한 기능과 연결된 미검증 문장이 없습니다.")
        st.warning("검증 절차: 해당 변형체 → 해당 실험군/비교군 → 효과 방향 → 정량값/p값 → 원저 식별정보의 일치 여부 확인 후 수동 또는 검증된 파이프라인으로 승격")

    if len(ranked) > 3:
        with st.expander(f"추가 검토 가능한 변형체 {len(ranked)-3}개 보기"):
            st.dataframe(pd.DataFrame([{"원형": show(x[3][0].get("Parent_Peptide")),
                                       "변형체": show(x[3][0].get("Variant_Name")),
                                       "독립 원저 식별": x[1]} for x in ranked[3:]]),
                         use_container_width=True, hide_index=True)

st.divider()
st.caption("문헌 원본: 기존 3차 DB + Europe PMC API. 새 논문은 build_database.py로 수집하며, Streamlit 화면은 SQLite를 읽기만 합니다. 동일 DOI 중복 결과와 검증 전 문장으로 추천 확신도를 높이지 않습니다.")
