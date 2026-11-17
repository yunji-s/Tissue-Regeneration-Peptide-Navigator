"""Peptide Modification Library v5 — one-click search, ingest, preprocess, display.

Run in the SAME private GitHub Codespace as the local Ollama server:
    python -m streamlit run streamlit_app_v5.py --server.address 0.0.0.0 --server.port 8502

Requires in this folder: navigator.sqlite3, local_evidence_labeler_v13.py.
Existing excel_screening is NEVER modified. New Europe PMC discoveries are
stored separately and DO NOT acquire inferred modification/sequence labels.
"""
from __future__ import annotations

import html
import json
import re
import sqlite3
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import streamlit as st

ROOT = Path(__file__).resolve().parent
DB = ROOT / 'navigator.sqlite3'
LABELER_FILE = ROOT / 'local_evidence_labeler_v13.py'
PEPTIDES = ('GHK', 'RADA16-I', 'QK', 'IKVAV', 'KLD-12', 'RGD', 'YIGSR')
API = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
PLACEHOLDERS = {'', 'nan', 'none', 'n/a', '-', '확인필요', '확인 필요', 'not reported'}
FUNCTIONS = {'WOUND':'상처 치유', 'ANGIO':'혈관신생', 'NERVE':'신경 재생',
             'CARTILAGE':'연골 재생', 'ADHESION':'세포 부착', 'MIGRATE':'세포 이동',
             'PROLIF':'세포 증식', 'ECMREMOD':'ECM 재형성', 'STABILITY':'안정성',
             'MYOGEN':'근육 형성', 'APOPTOSIS':'세포사멸', 'DIFF':'세포 분화',
             'DELIVERY':'전달', 'INFLAM':'염증', 'ANTIINFLAM':'항염증'}


def clean(value):
    return re.sub(r'\s+', ' ', str(value or '')).strip()


def show(value):
    s = clean(value)
    return '미확인' if s.lower() in PLACEHOLDERS else s


def doi_norm(value):
    return re.sub(r'^(?:https?://(?:dx\.)?doi\.org/|doi:\s*)', '', clean(value), flags=re.I).lower()


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def connect():
    conn = sqlite3.connect(DB, timeout=30)
    conn.row_factory = sqlite3.Row
    return conn


def ensure_discoveries(conn):
    # Only new, namespaced tables. Never rewrite Excel/Screening originals.
    conn.execute('''CREATE TABLE IF NOT EXISTS automatic_paper_discoveries_v1 (
       doi TEXT PRIMARY KEY, parent_peptide TEXT NOT NULL, title TEXT NOT NULL,
       journal TEXT NOT NULL, pub_year TEXT NOT NULL, abstract TEXT NOT NULL,
       paper_url TEXT NOT NULL, evidence_quote TEXT NOT NULL DEFAULT '',
       outcome_en TEXT NOT NULL DEFAULT '', evidence_scope TEXT NOT NULL DEFAULT '',
       evidence_status TEXT NOT NULL DEFAULT '', last_error TEXT NOT NULL DEFAULT '',
       processed_at TEXT NOT NULL DEFAULT '', source_hash TEXT NOT NULL DEFAULT ''
    )''')
    conn.commit()


def screening(conn):
    return [json.loads(row[0]) for row in conn.execute('SELECT payload_json FROM excel_screening')]


def existing_summaries(conn):
    tables = {r[0] for r in conn.execute('SELECT name FROM sqlite_master WHERE type="table"')}
    if 'qualitative_summaries_v2' not in tables:
        return {}
    return {r['record_id']: dict(r) for r in conn.execute('SELECT * FROM qualitative_summaries_v2')}


def ollama_available():
    try:
        with urlopen('http://127.0.0.1:11434/api/tags', timeout=2) as r:
            return r.status == 200
    except Exception:
        return False


def relevant_screening(rows, peptide, keyword):
    term = keyword.casefold().strip()
    hits = []
    for row in rows:
        original = clean(row.get('Original_Peptide'))
        if peptide != '7종 전체' and peptide not in [p.strip() for p in original.split('+')]:
            continue
        if term and term not in ' '.join(clean(v) for v in row.values()).casefold():
            continue
        hits.append(row)
    return hits


def epmc_query(peptide, keyword):
    # Quote peptide names, restrict to biomedical regeneration context; user keyword is a search term.
    base = '(regeneration OR wound OR cartilage OR nerve OR hydrogel OR angiogenesis OR tissue)'
    if peptide != '7종 전체':
        name = peptide.replace('"', '')
        base = f'"{name}" AND {base}'
    if keyword:
        safe = re.sub(r'["\[\]{}:()]', ' ', keyword).strip()
        if safe:
            base += f' AND ({" AND ".join(w for w in safe.split()[:5])})'
    return base


def epmc_search(peptide, keyword, limit=8):
    params = {'query': epmc_query(peptide, keyword), 'format':'json',
              'resultType':'core', 'pageSize':min(25, max(1,limit*3)), 'sort':'FIRST_PDATE_D desc'}
    req = Request(API + '?' + urlencode(params),
                  headers={'Accept':'application/json','User-Agent':'PeptideModificationLibrary/5.0 (research prototype)'})
    with urlopen(req, timeout=30) as resp:
        payload = json.load(resp)
    docs = []
    seen = set()
    for src in payload.get('resultList',{}).get('result',[]):
        doi = doi_norm(src.get('doi'))
        if not doi or doi in seen:
            continue
        title = clean(html.unescape(src.get('title') or ''))
        abstract = clean(html.unescape(re.sub(r'<[^>]+>', ' ', src.get('abstractText') or '')))
        if not abstract:
            continue
        # A search match is a candidate, not verification of modification-specific effect.
        seen.add(doi)
        docs.append({'doi':doi,'title':title,'abstract':abstract,
                     'journal':clean(src.get('journalTitle')),
                     'pub_year':clean(src.get('pubYear')),
                     'paper_url':'https://doi.org/'+doi})
        if len(docs) >= limit:
            break
    return docs


def short_summary(text):
    # Exactly up to 3 VISIBLE lines in the Streamlit layout.
    safe = html.escape(three_lines(text))
    st.markdown(
        '<div style="display:-webkit-box;-webkit-line-clamp:3;-webkit-box-orient:vertical;overflow:hidden;line-height:1.5;overflow-wrap:anywhere">'
        + safe + '</div>', unsafe_allow_html=True)


def three_lines(text):
    """Keep at most three substantive sentences, bounded text length."""
    s = clean(text)
    if not s:
        return '초록에서 보고된 결과를 확인하지 못했습니다.'
    pieces = re.split(r'(?<=[.!?])\s+(?=[A-Z(])', s)
    return ' '.join(pieces[:3])[:680]


def preprocess_new(conn, doc, peptide, labeler, use_llm):
    """Select an observed sentence, then optionally summarize. Do NOT invent sequence or modification."""
    doi = doc['doi']
    abstract = doc['abstract']
    import hashlib
    source_hash = hashlib.sha256((abstract+'|v5|'+peptide).encode('utf8')).hexdigest()
    old = conn.execute('SELECT * FROM automatic_paper_discoveries_v1 WHERE doi=?',(doi,)).fetchone()
    if old and old['source_hash'] == source_hash and old['processed_at']:
        return 'cached'
    # Use known original peptide token only as a relevance cue. No modification assignment.
    variant = peptide if peptide != '7종 전체' else ''
    ranked = labeler.rank_abstract_evidence(abstract, variant, metadata={})
    item = ranked[0] if ranked else None
    quote = item['quote'] if item else ''
    scope = item['scope'] if item else 'INSUFFICIENT_EVIDENCE'
    result = quote
    status = 'SOURCE_EXTRACT_REVIEW_REQUIRED' if quote else 'INSUFFICIENT_EVIDENCE'
    err = ''
    if quote and use_llm:
        try:
            candidate = labeler.summarize_one(labeler.MODEL, variant or 'studied materials', quote, scope,
                                               timeout=90, max_tokens=160)
            candidate, direction, issues = labeler.validate_summary(candidate, 'unclear', quote, variant or 'studied materials')
            if issues:
                err = 'Summary fidelity warnings: '+ '; '.join(issues)
            else:
                result = candidate
                status = 'LLM_AUTO_CHECKED_REVIEW_REQUIRED'
        except Exception as exc:
            err = f'{type(exc).__name__}: {str(exc)[:220]}'
    conn.execute('''INSERT INTO automatic_paper_discoveries_v1
      (doi,parent_peptide,title,journal,pub_year,abstract,paper_url,evidence_quote,outcome_en,
       evidence_scope,evidence_status,last_error,processed_at,source_hash)
      VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)
      ON CONFLICT(doi) DO UPDATE SET
       title=excluded.title,journal=excluded.journal,pub_year=excluded.pub_year,
       abstract=excluded.abstract,paper_url=excluded.paper_url,evidence_quote=excluded.evidence_quote,
       outcome_en=excluded.outcome_en,evidence_scope=excluded.evidence_scope,
       evidence_status=excluded.evidence_status,last_error=excluded.last_error,
       processed_at=excluded.processed_at,source_hash=excluded.source_hash''',
       (doi, peptide if peptide != '7종 전체' else '미분류',doc['title'],doc['journal'],doc['pub_year'],
        abstract,doc['paper_url'],quote,three_lines(result) if quote else '',scope,status,err,utc_now(),source_hash))
    conn.commit()
    return 'processed'


def preprocess_existing(record_ids, limit=8):
    """Invoke v13 automatically for matching rows only; never whole 57-record DB per search."""
    statuses = []
    for rid in list(dict.fromkeys(record_ids))[:limit]:
        if not LABELER_FILE.is_file():
            statuses.append((rid, 'v13 파일 없음'))
            continue
        try:
            res = subprocess.run([sys.executable,'-u',str(LABELER_FILE), '--db',str(DB),
                                  '--record-id',rid], capture_output=True,text=True,
                                 encoding='utf-8', errors='replace',timeout=150,cwd=ROOT)
            statuses.append((rid,'ok' if res.returncode==0 else (res.stderr or res.stdout)[-160:]))
        except subprocess.TimeoutExpired:
            statuses.append((rid,'timeout'))
    return statuses


def display_existing(row, summary):
    name = show(row.get('Modified_Peptide'))
    with st.container(border=True):
        st.markdown('#### '+name)
        c1,c2 = st.columns(2)
        with c1:
            st.write('**원형 펩타이드:**', show(row.get('Original_Peptide')))
            st.write('**원형 서열:**', show(row.get('Original_Sequence')))
            st.write('**Modification 종류:**', show(row.get('Modification_Method')))
            st.write('**Modification 위치:**', show(row.get('Modification_Position')))
        with c2:
            st.write('**변형 서열/구조:**', show(row.get('Modified_Sequence')))
            st.write('**Formulation:**', show(row.get('Formulation')))
            st.write('**연구 기능:**',show(row.get('Target_Function')))
            st.write('**평가 지표:**',show(row.get('Measured_Outcome')))
        st.markdown('**논문 요약 (최대 3문장)**')
        if summary and clean(summary.get('summary_ko')):
            short_summary(summary['summary_ko'])
            st.caption('AI 자동 요약 / 초록 근거 · 사람의 검토 전 · '+show(summary.get('verification_status')))
        else:
            st.write('검증 가능한 요약이 아직 없습니다. 기존 DB 원문 정보를 유지합니다.')
        st.caption('Record ID: '+show(row.get('Record_ID'))+' · '+show(row.get('Reference')))
        doi=doi_norm(row.get('DOI'))
        url=clean(row.get('Paper_URL'))
        if not url.startswith(('https://','http://')):
            url='https://doi.org/'+doi if doi else ''
        if url:
            st.link_button('연결 논문 열기',url)
        else:
            st.caption('논문 링크 미확인')


def display_new(doc):
    with st.container(border=True):
        st.markdown('#### '+show(doc['title']))
        c1,c2=st.columns(2)
        with c1:
            st.write('**검색 기준 펩타이드:**',show(doc['parent_peptide']))
            st.write('**Modification 종류:**','미확인 (원문 검토 필요)')
            st.write('**Modification 위치:**','미확인')
        with c2:
            st.write('**변형 서열/구조:**','미확인')
            st.write('**Formulation:**','미확인')
            st.write('**연구 기능:**','미확인 (초록만으로 확정 불가)')
        st.markdown('**논문 요약 (최대 3문장)**')
        short_summary(doc['outcome_en'] if doc['outcome_en'] else '초록에서 결과 근거를 확인하지 못했습니다.')
        st.caption('신규 수집 · 자동 분석 / 미검증 · '+show(doc['evidence_status']))
        if doc['last_error']:
            with st.expander('분석 검토 메모'):
                st.write(doc['last_error'])
        st.link_button('연결 논문 열기',doc['paper_url'])


def execute_search(peptide, term, include_new, ollama, max_new=6):
    with connect() as conn:
        ensure_discoveries(conn)
        rows = relevant_screening(screening(conn),peptide,term)
    unique_ids = [r['Record_ID'] for r in rows if doi_norm(r.get('DOI'))]
    # Preprocessing is initiated by this SEARCH, not a separate button.
    logs=[]
    if ollama and unique_ids:
        logs=preprocess_existing(unique_ids,limit=8)
    elif unique_ids and not ollama:
        logs=[('Ollama','연결 불가: 기존 DB 요약을 사용')]
    new_dois=[]
    warning=''
    if include_new:
        try:
            papers=epmc_search(peptide,term,limit=max_new)
            with connect() as conn:
                ensure_discoveries(conn)
                known_dois={doi_norm(r.get('DOI')) for r in screening(conn) if doi_norm(r.get('DOI'))}
                for paper in papers:
                    if paper['doi'] in known_dois:
                        continue
                    preprocess_new(conn,paper,peptide,sys.modules['local_evidence_labeler_v13'],ollama)
                    new_dois.append(paper['doi'])
        except Exception as exc:
            warning=f'Europe PMC 검색 또는 신규 논문 분석 중 오류: {type(exc).__name__}: {exc}'
    return {'rows_ids':[r['Record_ID'] for r in rows], 'new_dois':new_dois, 'logs':logs,
            'warning':warning,'peptide':peptide,'term':term,'done_at':utc_now()}


st.set_page_config(page_title='Peptide Modification Library',page_icon='🧬',layout='wide')
st.title('🧬 Peptide Modification Library')
st.caption('검색과 동시에 기존 DB 연결·논문 수집·결과 전처리를 수행합니다. 원본 Screening DB 65건은 보존됩니다.')
if not DB.is_file():
    st.error('navigator.sqlite3가 없습니다. Codespaces 프로젝트 폴더에 최신 DB를 넣어주십시오.')
    st.stop()
if not LABELER_FILE.is_file():
    st.error('local_evidence_labeler_v13.py가 같은 폴더에 없습니다.')
    st.stop()
# Import existing deterministic evidence selector and validation logic, rather than copying it.
import local_evidence_labeler_v13  # noqa: E402

with st.form('research_search'):
    c1,c2=st.columns([1,2])
    with c1:
        peptide=st.selectbox('원형 펩타이드',('7종 전체',)+PEPTIDES,index=0)
    with c2:
        term=st.text_input('검색어 (변형체·연구 기능·소재·논문 키워드)',
                           placeholder='예: GHK-Cu / hydrogel / nerve regeneration')
    include_new=st.checkbox('Europe PMC에서 신규 논문도 수집',value=True)
    submitted=st.form_submit_button('🔍 검색 및 자동 전처리',type='primary',use_container_width=True)

if submitted:
    is_ollama=ollama_available()
    with st.spinner('기존 기록 조회 → 신규 논문 확보 → 근거 추출 및 요약을 처리 중입니다…'):
        result=execute_search(peptide,term.strip(),include_new,is_ollama)
        st.session_state['latest_search']=result
    if result['warning']:
        st.warning(result['warning'])
    if not is_ollama:
        st.info('Ollama 미연결: 기존 DB 요약 또는 논문 원문 결과 문장을 사용했습니다. 신규 논문에 AI 요약을 생성하지 않았습니다.')
    if len(result['rows_ids'])>8 and is_ollama:
        st.info('화면 검색 결과는 모두 표시하지만, 요청당 기존 기록 전처리는 최대 8건입니다. 나머지는 저장된 요약을 표시합니다.')

if 'latest_search' not in st.session_state:
    st.info('검색어를 입력한 후 검색 버튼을 누르면 자동으로 전처리하고 결과를 표시합니다. 검색어 없이도 기존 7종 DB를 조회할 수 있습니다.')
else:
    result=st.session_state['latest_search']
    with connect() as conn:
        ensure_discoveries(conn)
        original={r.get('Record_ID'):r for r in screening(conn)}
        summaries=existing_summaries(conn)
        new=[]
        for doi in result['new_dois']:
            row=conn.execute('SELECT * FROM automatic_paper_discoveries_v1 WHERE doi=?',(doi,)).fetchone()
            if row:
                new.append(dict(row))
    hits=[original[rid] for rid in result['rows_ids'] if rid in original]
    st.subheader('검색 결과')
    m1,m2,m3=st.columns(3)
    m1.metric('기존 Screening 기록',len(hits))
    m2.metric('새 논문 (중복 제외)',len(new))
    m3.metric('연결 논문 DOI',len({doi_norm(r.get('DOI')) for r in hits if doi_norm(r.get('DOI'))} | {r['doi'] for r in new}))
    st.caption('기존 DB 필드는 그대로 표시합니다. 신규 논문의 Modification 종류·서열·기능은 근거 없이 자동 확정하지 않습니다.')
    if result['warning']:
        st.warning(result['warning'])
    if hits:
        st.markdown('### 기존 DB의 Modification')
        for row in hits:
            display_existing(row,summaries.get(row['Record_ID']))
    if new:
        st.markdown('### 새로 수집한 논문 (Modification 검토 필요)')
        for row in new:
            display_new(row)
    if not hits and not new:
        st.info('일치하는 기록을 찾지 못했습니다. 검색어를 바꾸거나 원형 펩타이드를 선택해 보십시오.')
    if result['logs']:
        with st.expander('자동 전처리 실행 정보'):
            for rid,status in result['logs']:
                st.write(f'{rid}: {status}')
    st.caption('새로 수집된 논문은 별도 테이블에 저장됩니다. Codespaces의 SQLite 변경 사항은 GitHub에 자동 커밋되지 않습니다.')
