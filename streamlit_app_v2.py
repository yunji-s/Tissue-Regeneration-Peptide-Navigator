"""Read-only Tissue Regeneration Peptide Navigator UI.
Run data building separately: python build_database.py --peptide GHK
Launch: streamlit run streamlit_app_v2.py
"""
import json
import sqlite3
from pathlib import Path
import streamlit as st
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATABASE = ROOT / 'navigator.sqlite3'
FUNCTIONS = {'WOUND':'Wound healing', 'ANGIO':'Angiogenesis','NERVE':'Nerve regeneration'}
PEPTIDES = {'GHK':'P001','RADA16-I':'P002','QK':'P003','IKVAV':'P004','KLD-12':'P005','RGD':'P006','YIGSR':'P007'}

def valid(s):
    return str(s or '').strip()

def source_url(row):
    url=valid(row.get('Paper_URL'))
    if url.startswith(('https://','http://')):return url
    doi=valid(row.get('DOI'))
    return 'https://doi.org/'+doi if doi else ''

def load_data(peptide,code):
    with sqlite3.connect(DATABASE) as db:
        recs = [json.loads(x[0]) for x in db.execute('''SELECT payload_json FROM excel_recommendations WHERE peptide_id=? AND function_code=?''',(PEPTIDES[peptide],code))]
        old = [json.loads(x[0]) for x in db.execute('''SELECT payload_json FROM excel_screening WHERE peptide=? AND function_code=?''',(peptide,code))]
        new = [json.loads(x[0]) for x in db.execute('''SELECT e.screening_json FROM extraction_candidates e WHERE e.peptide=?''',(peptide,)) if json.loads(x[0]).get('Target_Function_Code') == code]
        stats = {t:db.execute('SELECT count(*) FROM '+t).fetchone()[0] for t in ('papers','extraction_candidates')}
    return recs,old,new,stats

def tier(row):
    status=valid(row.get('Evidence_Status')).lower()
    eligibility=valid(row.get('Recommendation_Eligibility')).lower()
    if 'unverified' in status or 'ineligible' in eligibility or 'context only' in eligibility:
        return '탐색/맥락 근거'
    if 'pending' in status or 'pending' in eligibility or 'caveat' in eligibility:
        return '보고된 결과 · 원문 검증 필요'
    # The historic Excel status alone is NOT sufficient to confirm an actual audit.
    return '기존 DB 기록 · 원문 검증 상태 재확인 필요'

def rank_key(row):
    # Conservative DISPLAY ORDER, not a scientific efficacy score.
    category=tier(row)
    priority=0 if category.startswith('기존 DB') else 1 if category.startswith('보고된') else 2
    comparator=valid(row.get('Comparator_Type')).lower()
    comparator_priority=0 if any(x in comparator for x in ('direct','parent','original')) else 1
    return (priority,comparator_priority,valid(row.get('Recommendation_ID')))

st.set_page_config(page_title='Tissue Regeneration Peptide Navigator',page_icon='🧬',layout='wide')
st.title('🧬 Tissue Regeneration Peptide Navigator')
st.caption('기존 3차 DB 기반 기능별 문헌 후보 추천 · Europe PMC 자동 추출은 별도 구축 프로세스에서 수행')
if not DATABASE.is_file():
    st.error('navigator.sqlite3가 없습니다. 먼저 `python build_database.py --report`를 실행하세요.')
    st.stop()
left,right = st.columns(2)
with left: peptide=st.selectbox('펩타이드',list(PEPTIDES))
with right: function=st.selectbox('목표 기능',list(FUNCTIONS),format_func=lambda x:FUNCTIONS[x])
if st.button('기존 문헌 근거 조회',type='primary',use_container_width=True):
    try: recs,old,new,stats=load_data(peptide,function)
    except (sqlite3.DatabaseError,KeyError,ValueError) as e:
        st.error(f'DB 조회 실패: {e}');st.stop()
    st.info(f'원본 Screening 기록 {len(old)}건 · 기존 추천 기록 {len(recs)}건 · 새 자동추출 문장 {len(new)}건 (검증 전)')
    st.subheader('기능별 후보 · 최대 3개')
    unique={}
    for row in sorted(recs,key=rank_key):
        k=(valid(row.get('Variant_Name')),valid(row.get('Modification_Method')),valid(row.get('Formulation_Method')))
        if k not in unique:unique[k]=row
    eligible=[r for r in unique.values() if tier(r)!='탐색/맥락 근거'][:3]
    if not eligible: st.warning('이 기능에서 표시 가능한 기존 추천 근거가 없습니다. 검증되지 않은 후보로 순위를 채우지 않습니다.')
    for i,r in enumerate(eligible,1):
        with st.container(border=True):
            st.markdown(f"**{i}순위 (잠정) — {valid(r.get('Variant_Name')) or peptide}**")
            st.write('**근거 상태:**',tier(r))
            st.write('**변형:**',valid(r.get('Modification_Method')) or 'Not reported')
            st.write('**제형:**',valid(r.get('Formulation_Method')) or 'Not reported')
            st.write('**평가 지표:**',valid(r.get('Outcome_Metric')) or 'Not reported')
            st.write('**효과 방향(기존 DB):**',valid(r.get('Observed_Direction')) or 'Not reported')
            st.write('**비교군:**',valid(r.get('Comparator')) or 'Not reported')
            st.write('**한계:**',valid(r.get('Limitations')) or '원문 확인 필요')
            if source_url(r):st.markdown(f'[원문/DOI]({source_url(r)})')
    st.caption('순위는 검증 상태와 비교군 유형을 기준으로 한 임시 표시 순서이며, 실제 효능 우위 또는 통계적으로 검증된 추천 점수가 아닙니다.')
    with st.expander('탐색/맥락 근거 및 기존 Screening 원자료'):
        if old: st.dataframe(pd.DataFrame(old),use_container_width=True,hide_index=True)
        else: st.write('해당 기능의 기존 Screening 기록 없음')
        excluded=[r for r in recs if tier(r)=='탐색/맥락 근거']
        if excluded:st.write(f'추천 미반영 맥락 근거 {len(excluded)}건')
    st.subheader('새로 자동 수집한 미검증 문장')
    st.warning('아래 문장의 증가/감소, 정량값, 실험군·비교군은 검증되지 않았습니다. 추천에 반영하지 않습니다.')
    if new:
        cols=['Original_Peptide','Target_Function_Code','Candidate_Statement','Source_Verification','DOI','Evidence_Status']
        st.dataframe(pd.DataFrame(new).reindex(columns=cols),use_container_width=True,hide_index=True)
    else:st.info('아직 수집된 후보가 없습니다. Europe PMC 수집은 build_database.py에서 별도로 수행합니다.')
