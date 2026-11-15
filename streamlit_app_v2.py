"""Minimal Streamlit recommendation display; no AI API and no network calls."""
from pathlib import Path
import streamlit as st
from recommendation_engine import FUNCTIONS, PEPTIDES, load, recommend

ROOT = Path(__file__).resolve().parent
DB = ROOT / 'navigator.sqlite3'
st.set_page_config(page_title='Tissue Regeneration Peptide Navigator', page_icon='🧬', layout='wide')
st.title('🧬 Tissue Regeneration Peptide Navigator')
st.caption('7종 펩타이드 기반 Modification·Formulation 근거 검토 프로토타입')

if not DB.is_file():
    st.error('navigator.sqlite3 파일이 없습니다. Codespaces에서 구축 후 GitHub에 커밋하십시오.')
    st.stop()

@st.cache_data(show_spinner=False)
def get_data(path, modified, size):
    return load(path)

try:
    stat = DB.stat()
    recommendations, screening, processed = get_data(str(DB), stat.st_mtime_ns, stat.st_size)
except Exception as exc:
    st.error(f'데이터베이스 조회 오류: {exc}')
    st.stop()

with st.form('selector'):
    a, b = st.columns([3, 2])
    with a:
        function = st.selectbox('목표 조직재생 기능', list(FUNCTIONS), format_func=lambda k: FUNCTIONS[k])
    with b:
        scope = st.selectbox('검색 범위', ['7종 전체', *PEPTIDES])
    submitted = st.form_submit_button('Modification·Formulation 후보 검색', type='primary', use_container_width=True)

if submitted:
    choices = recommend(recommendations, screening, processed, function, scope)
    st.subheader('추천 검토 후보')
    st.caption('논문의 정량 효과를 서로 직접 비교한 순위가 아니라 기존 DB의 근거 상태·연결 DOI 수에 따른 검토 우선순위입니다.')
    if not choices:
        st.info('해당 조건에서 확인 가능한 기존 DB 근거가 부족합니다. 후보를 임의로 추가하지 않습니다.')
    for rank, entry in enumerate(choices, 1):
        with st.container(border=True):
            st.markdown(f"### {rank}. {entry['variant']}  ·  원형 {entry['parent']}")
            st.write('**추천 Modification:**', entry['modification'])
            st.write('**Formulation:**', entry['formulation'])
            st.write('**변형 서열/구조:**', entry['sequence'])
            st.write('**문헌 기반 기대 기능:**', entry['function'])
            st.write('**근거 요약:**', entry['summary'])
            st.caption(f"연결 DOI {entry['doi_count']}개 · 독립 연구 확인 필요 · 자동 전처리 문장 {entry['processed_mentions']}건 연결(추천 점수 미반영)")
            # Keep links short; no raw evidence dumps or long preprocessing tables.
            if entry['links']:
                first = next((e for e in entry['links'] if e['url']), None)
                if first:
                    st.markdown(f"[대표 원저 열기]({first['url']})")
                with st.expander('관련 원저 링크'):
                    seen = set()
                    for paper in entry['links']:
                        if paper['doi'] in seen or not paper['url']:
                            continue
                        seen.add(paper['doi'])
                        st.markdown(f"[{paper['reference'] or paper['doi']}]({paper['url']})")

st.divider()
st.caption('AI API 미사용 · DB 기반 규칙형 요약 · 원문 문장/정량값은 전처리 DB에 보존되며 검증 없이 효능 사실로 확정하지 않습니다.')
