"""Evidence-aware 7-peptide recommendation helper (no AI API).

Read-only SQLite. Literature-derived candidate priority, NOT validated efficacy ranking.
Processed evidence is attached transparently, never used to raise recommendation rank.
"""
from __future__ import annotations
import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

PEPTIDES = ('GHK','RADA16-I','QK','IKVAV','KLD-12','RGD','YIGSR')
FUNCTIONS = {'WOUND':'Wound healing','ANGIO':'Angiogenesis','NERVE':'Nerve regeneration',
             'CARTILAGE':'Cartilage regeneration','CART':'Cartilage regeneration',
             'INFLAM':'Inflammatory response','ANTIINFLAM':'Anti-inflammation',
             'ADHESION':'Cell adhesion','MIGRATE':'Cell migration',
             'PROLIF':'Cell proliferation','ECMREMOD':'ECM remodeling'}
BAD = ('unverified','context only','ineligible','exclude','recheck')
NEG = ('non-significant','nonsignificant','not significant','no significant')
NULLS = ('','none','nan','not reported','확인필요','확인 필요','-')


def value(x): return str(x or '').strip()
def doi(x): return re.sub(r'^https?://(?:dx\.)?doi\.org/','',value(x),flags=re.I).lower()
def missing(x): return value(x).lower() in NULLS


def load(path):
    path = Path(path)
    if not path.is_file(): raise FileNotFoundError(f'SQLite 파일이 없습니다: {path}')
    with sqlite3.connect(f'file:{path.resolve()}?mode=ro',uri=True) as db:
        tables = {x[0] for x in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        required = {'excel_recommendations','excel_screening'}
        if required - tables: raise ValueError('누락된 DB 테이블: '+', '.join(sorted(required-tables)))
        rec = [json.loads(x[0]) for x in db.execute('SELECT payload_json FROM excel_recommendations')]
        screening = [json.loads(x[0]) for x in db.execute('SELECT payload_json FROM excel_screening')]
        processed = []
        if 'processed_evidence' in tables:
            for key,pep,fun,d,payload,status in db.execute('''SELECT paper_key,peptide,function_code,doi,
                         result_json,extraction_status FROM processed_evidence'''):
                try: item=json.loads(payload or '{}')
                except (TypeError,ValueError): continue
                if not isinstance(item,dict): continue
                item.update(paper_key=key,peptide=pep,function_code=fun,doi=doi(d),extraction_status=status)
                processed.append(item)
    return rec,screening,processed


def eligible(row):
    state = ' | '.join(value(row.get(f)).lower() for f in ('Evidence_Status','Recommendation_Eligibility','Observed_Direction'))
    return not any(word in state for word in BAD+NEG)


def result_note(rows):
    # Preserve the source's wording rather than generating a causal efficacy claim.
    for row in rows:
        content = value(row.get('Evidence_Detail'))
        if not missing(content): return content
    for row in rows:
        content = value(row.get('Short_Summary'))
        if not missing(content): return content
    return '기존 DB에 기능 관련 기록이 있으나 구체적인 결과 설명이 없습니다.'


def record_url(row):
    u=value(row.get('Paper_URL'))
    d=doi(row.get('DOI'))
    return u if u.startswith(('http://','https://')) else ('https://doi.org/'+d if d else '')


def recommend(recommendations,screening,processed,function,scope='7종 전체',limit=3):
    wanted = set(PEPTIDES) if scope=='7종 전체' else {scope}
    grouped=defaultdict(list)
    for r in recommendations:
        if value(r.get('Parent_Peptide')) not in wanted or value(r.get('Target_Function_Code'))!=function: continue
        key=tuple(value(r.get(f)).casefold() for f in ('Parent_Peptide','Variant_Name','Modification_Method','Formulation_Method'))
        grouped[key].append(r)
    results=[]
    for key,rows in grouped.items():
        usable=[r for r in rows if eligible(r)]
        if not usable: continue
        base=usable[0]
        parent=value(base.get('Parent_Peptide'))
        variant=value(base.get('Variant_Name'))
        ids={value(r.get('Variant_ID')) for r in usable if value(r.get('Variant_ID'))}
        dois={doi(r.get('DOI')) for r in usable if doi(r.get('DOI'))}
        related=[s for s in screening if doi(s.get('DOI')) in dois and
                 ((value(s.get('Variant_ID')) in ids and bool(ids)) or
                  value(s.get('Modified_Peptide')).casefold()==variant.casefold())]
        seq=next((value(s.get('Modified_Sequence')) for s in related if not missing(s.get('Modified_Sequence'))),'Not reported')
        evidence=[]
        for r in usable:
            desc=result_note([r])
            evidence.append({'type':'기존 3차 DB','doi':doi(r.get('DOI')),'url':record_url(r),
              'title':value(r.get('Reference')),'text':desc,
              'comparator':value(r.get('Comparator')),'metric':value(r.get('Outcome_Metric')),
              'status':value(r.get('Evidence_Status'))})
        direct=[]
        indirect=[]
        for p in processed:
            if doi(p.get('doi')) not in dois or value(p.get('peptide')).casefold()!=parent.casefold(): continue
            if value(p.get('function_code')) and value(p.get('function_code'))!=function: continue
            mentioned=value(p.get('matched_variant'))
            if missing(mentioned): continue
            sentence=value(p.get('source_sentence'))
            if not sentence: continue
            entry={'type':'Europe PMC 전처리 · 미검증','doi':doi(p.get('doi')),
                   'url':'https://doi.org/'+doi(p.get('doi')) if doi(p.get('doi')) else '',
                   'title':'논문 추출 문장','text':sentence,
                   'status':value(p.get('extraction_status')),'metric':value(p.get('metric')),
                   'relative':p.get('relative_change_percent'),
                   'reported_relative':p.get('reported_relative_change_percent'),
                   'review_reason':value(p.get('review_reason'))}
            if mentioned.casefold()==variant.casefold(): direct.append(entry)
            elif mentioned.casefold()==parent.casefold(): indirect.append(entry)
        # Deduplicate same sentences/rows, preserving provenance
        def unique(entries):
            seen=set(); out=[]
            for e in entries:
                k=(e['type'],e['doi'],e['text'])
                if k not in seen: seen.add(k); out.append(e)
            return out
        evidence=unique(evidence)
        direct=unique(direct)
        indirect=unique(indirect)
        marked=any('verified significant' in value(r.get('Evidence_Status')).lower() and
                   'eligible (verified)' in value(r.get('Recommendation_Eligibility')).lower() for r in usable)
        note=result_note(usable)
        brief=(f'기존 DB에는 “{note[:135]}{"…" if len(note)>135 else ""}”라고 기록됨. '
               f'연결 DOI {len(dois)}개이며, 변형체 직접 비교 및 독립 재현은 추가 확인이 필요합니다.')
        results.append({'parent':parent,'variant':variant,'modification':value(base.get('Modification_Method')) or 'Not reported',
           'formulation':value(base.get('Formulation_Method')) or 'Not reported','sequence':seq,
           'function':FUNCTIONS.get(function,function),'summary':brief,'doi_count':len(dois),
           'marked_verified_in_excel':marked,'evidence':evidence,'direct_excerpts':direct,
           'parent_context_excerpts':indirect,'all_preprocessed_matches':len(direct)+len(indirect)})
    results.sort(key=lambda r:(-int(r['marked_verified_in_excel']),-r['doi_count'],r['parent'].lower(),r['variant'].lower()))
    return results[:limit]


if __name__=='__main__':
    r,s,p=load(Path(__file__).with_name('navigator.sqlite3'))
    print(f'Loaded: recommendation {len(r)}, screening {len(s)}, processed {len(p)}')
    for f in ('WOUND','ANGIO','NERVE'):
        print(f,[(x['parent'],x['variant'],len(x['evidence']),len(x['direct_excerpts'])) for x in recommend(r,s,p,f)])
