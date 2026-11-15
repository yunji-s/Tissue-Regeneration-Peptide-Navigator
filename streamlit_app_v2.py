"""Concise recommendation cards + optional cited evidence (offline, no LLM API)."""
from pathlib import Path
import streamlit as st
from recommendation_engine import FUNCTIONS,PEPTIDES,load,recommend

ROOT=Path(__file__).resolve().parent
DB=ROOT/'navigator.sqlite3'
st.set_page_config(page_title='Tissue Regeneration Peptide Navigator',page_icon='🧬',layout='wide')
st.title('🧬 Tissue Regeneration Peptide Navigator')
st.caption('7종 펩타이드 · Modification·Formulation 근거 검토 · AI API 미사용')

if not DB.is_file():
    st.error('navigator.sqlite3 파일이 없습니다. Codespaces의 DB를 GitHub에 Push하십시오.')
    st.stop()

@st.cache_data(show_spinner=False)
def read_db(path,mtime,size): return load(path)

try:
    stat=DB.stat()
    recommendations,screenings,processed=read_db(str(DB),stat.st_mtime_ns,stat.st_size)
except Exception as exc:
    st.error('DB 조회 오류: '+str(exc))
    st.stop()

with st.form('filters'):
    left,right=st.columns([3,2])
    with left:
        function=st.selectbox('목표 조직재생 기능',list(FUNCTIONS),format_func=lambda k:FUNCTIONS[k])
    with right:
        scope=st.selectbox('검색 범위',['7종 전체',*PEPTIDES])
    submitted=st.form_submit_button('근거 기반 후보 검색',type='primary',use_container_width=True)

if 'active_search' not in st.session_state: st.session_state.active_search=None
if submitted: st.session_state.active_search=(function,scope)

if st.session_state.active_search:
    f,s=st.session_state.active_search
    choices=recommend(recommendations,screenings,processed,f,s)
    st.subheader('Modification·Formulation 검토 후보')
    st.caption('검토 우선순위이며 효능 우열이 검증된 순위가 아닙니다. 신규 자동 추출 결과는 추천 순위를 변경하지 않습니다.')
    if not choices: st.info('해당 기능의 추천 근거가 충분하지 않습니다.')
    for i,c in enumerate(choices,1):
        with st.container(border=True):
            st.markdown(f"### {i}. {c['variant']} · {c['parent']}")
            st.write('**추천 Modification:**',c['modification'])
            st.write('**Formulation:**',c['formulation'])
            st.write('**변형 서열/구조:**',c['sequence'])
            st.write('**문헌 기반 기대 기능:**',c['function'])
            st.write('**근거 요약:**',c['summary'])
            st.caption(f"기존 DB 연결 DOI {c['doi_count']}개 · 전처리된 직접 명칭 일치 문장 {len(c['direct_excerpts'])}건 (미검증)")
            with st.expander('근거 더보기 · 원문과 수치 확인'):
                if c['evidence']:
                    st.markdown('**기존 DB에 기록된 근거**')
                    for j,e in enumerate(c['evidence'][:5],1):
                        st.write(f"{j}. {e['text']}")
                        if e['metric']: st.caption('평가 지표: '+e['metric'])
                        if e['comparator']: st.caption('비교군: '+e['comparator'])
                        if e['url']: st.markdown(f"[근거 원저 {j} 열기]({e['url']})")
                if c['direct_excerpts']:
                    st.markdown('**Europe PMC 자동 추출 문장 — 변형체 명칭 직접 일치 (검증 전)**')
                    for j,e in enumerate(c['direct_excerpts'][:4],1):
                        st.write(f'원문에 “{e["text"]}”라고 명시됨')
                        if e['relative'] is not None:
                            st.caption(f"자동 계산 상대 변화율 {e['relative']:+g}% (원저 확인 전)")
                        elif e['reported_relative'] is not None:
                            st.caption(f"원문 보고 상대 변화율 {e['reported_relative']:+g}% (원저 확인 전)")
                        if e['url']: st.markdown(f"[추출 원저 {j} 열기]({e['url']})")
                if c['parent_context_excerpts']:
                    st.markdown('**원형 펩타이드 명칭만 확인된 문장 (변형체 직접 근거 아님)**')
                    st.caption(f"{len(c['parent_context_excerpts'])}건: 원형과 변형체의 효과를 혼동하지 않도록 기본 화면에서 제외합니다.")
                if not c['evidence'] and not c['direct_excerpts']:
                    st.info('연결된 문헌 근거가 없습니다.')

st.divider()
st.caption(f'현재 로드한 자동 전처리 결과 {len(processed)}건 · 직접 연결되지 않은 문장도 SQLite에 보존됩니다. '
           '전처리는 원저 검증을 대신하지 않으며, 실험 조건이 다른 결과를 직접 비교하지 않습니다.')
