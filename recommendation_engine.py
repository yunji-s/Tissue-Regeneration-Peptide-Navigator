"""Evidence-first, offline recommendation helper; no AI API or network calls.

Reads existing navigator.sqlite3. No changes to Excel/SQLite.
This is a screening-priority list, not a validated efficacy ranking.
"""
import json
import re
import sqlite3
from collections import defaultdict
from pathlib import Path

PEPTIDES = ('GHK', 'RADA16-I', 'QK', 'IKVAV', 'KLD-12', 'RGD', 'YIGSR')
FUNCTIONS = {'WOUND':'Wound healing', 'ANGIO':'Angiogenesis', 'NERVE':'Nerve regeneration',
             'CARTILAGE':'Cartilage regeneration', 'CART':'Cartilage regeneration',
             'INFLAM':'Inflammatory response', 'ANTIINFLAM':'Anti-inflammation',
             'ADHESION':'Cell adhesion', 'MIGRATE':'Cell migration',
             'PROLIF':'Cell proliferation', 'ECMREMOD':'ECM remodeling'}
BAD = ('unverified', 'context only', 'ineligible', 'exclude', 'recheck')
NEG = ('non-significant', 'nonsignificant', 'not significant', 'no significant')


def value(item):
    return str(item or '').strip()


def doi(item):
    return re.sub(r'^https?://(?:dx\.)?doi\.org/', '', value(item), flags=re.I).lower()


def missing(item):
    return value(item).lower() in ('', 'none', 'nan', 'not reported', '확인필요', '확인 필요', '-')


def load(path):
    if not Path(path).is_file():
        raise FileNotFoundError(f'SQLite 파일이 없습니다: {path}')
    with sqlite3.connect(f'file:{Path(path).resolve()}?mode=ro', uri=True) as conn:
        tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        required = {'excel_recommendations', 'excel_screening'}
        if required - tables:
            raise ValueError('누락된 DB 테이블: ' + ', '.join(sorted(required - tables)))
        recommendations = [json.loads(r[0]) for r in conn.execute('SELECT payload_json FROM excel_recommendations')]
        screening = [json.loads(r[0]) for r in conn.execute('SELECT payload_json FROM excel_screening')]
        processed = []
        if 'processed_evidence' in tables:
            for paper_key, peptide, function_code, d, payload in conn.execute(
                'SELECT paper_key, peptide, function_code, doi, result_json FROM processed_evidence'):
                try:
                    item = json.loads(payload)
                except (TypeError, ValueError):
                    continue
                item.update(paper_key=paper_key, peptide=peptide, function_code=function_code, doi=doi(d))
                processed.append(item)
    return recommendations, screening, processed


def eligibility(row):
    status = value(row.get('Evidence_Status')).lower()
    elig = value(row.get('Recommendation_Eligibility')).lower()
    outcome = value(row.get('Observed_Direction')).lower()
    if any(x in status or x in elig or x in outcome for x in NEG + BAD):
        return False
    return True


def summary(rows, source_count):
    metrics = sorted({value(r.get('Outcome_Metric')) for r in rows if not missing(r.get('Outcome_Metric'))})
    directions = sorted({value(r.get('Observed_Direction')) for r in rows if not missing(r.get('Observed_Direction'))})
    first_metric = ', '.join(metrics[:2]) if metrics else '평가 지표'
    first_direction = ', '.join(directions[:2]) if directions else '효과 방향 미확인'
    if source_count == 1:
        quantity = '연결된 원저 DOI는 1개이며'
    elif source_count:
        quantity = f'연결된 서로 다른 DOI는 {source_count}개이며'
    else:
        quantity = '확인된 DOI가 없으며'
    return (f'기존 DB에 {first_metric} 지표에서 {first_direction} 결과가 기록되어 있습니다. '
            f'{quantity}, 실험 조건·변형체 직접 비교·독립 재현 여부는 추가 확인이 필요합니다.')


def recommend(recommendations, screening, processed, function, scope='7종 전체', limit=3):
    names = set(PEPTIDES) if scope == '7종 전체' else {scope}
    groups = defaultdict(list)
    for r in recommendations:
        if value(r.get('Parent_Peptide')) not in names or value(r.get('Target_Function_Code')) != function:
            continue
        key = tuple(value(r.get(k)).casefold() for k in
                    ('Parent_Peptide', 'Variant_Name', 'Modification_Method', 'Formulation_Method'))
        groups[key].append(r)
    candidates = []
    for key, all_rows in groups.items():
        usable = [r for r in all_rows if eligibility(r)]
        if not usable:
            continue
        rep = usable[0]
        related_dois = {doi(r.get('DOI')) for r in usable if doi(r.get('DOI'))}
        matched_screen = [s for s in screening if doi(s.get('DOI')) in related_dois and
                          value(s.get('Variant_ID')) in {value(r.get('Variant_ID')) for r in usable}]
        seq = next((value(s.get('Modified_Sequence')) for s in matched_screen
                    if not missing(s.get('Modified_Sequence'))), 'Not reported')
        # Link only on DOI plus exact variant/parent. Never use unmatched snippets as proof.
        variant = value(rep.get('Variant_Name')).casefold()
        parent = value(rep.get('Parent_Peptide')).casefold()
        linked = [p for p in processed if p.get('doi') in related_dois and
                  value(p.get('peptide')).casefold() == parent and
                  value(p.get('matched_variant')).casefold() == variant]
        marked = any('verified significant' in value(r.get('Evidence_Status')).lower()
                     and 'eligible (verified)' in value(r.get('Recommendation_Eligibility')).lower()
                     for r in usable)
        candidates.append({
            'parent': value(rep.get('Parent_Peptide')),
            'variant': value(rep.get('Variant_Name')),
            'modification': value(rep.get('Modification_Method')) or 'Not reported',
            'formulation': value(rep.get('Formulation_Method')) or 'Not reported',
            'sequence': seq,
            'function': FUNCTIONS.get(function, function),
            'summary': summary(usable, len(related_dois)),
            'doi_count': len(related_dois),
            'marked_verified_in_excel': marked,
            'links': [{'doi': doi(r.get('DOI')), 'url': value(r.get('Paper_URL')) or
                       ('https://doi.org/'+doi(r.get('DOI')) if doi(r.get('DOI')) else ''),
                       'reference':value(r.get('Reference')), 'metric':value(r.get('Outcome_Metric'))}
                      for r in usable],
            'processed_mentions': len(linked),
        })
    candidates.sort(key=lambda r: (-int(r['marked_verified_in_excel']), -r['doi_count'],
                                   r['parent'].lower(), r['variant'].lower()))
    return candidates[:limit]


if __name__ == '__main__':
    r, s, p = load(Path(__file__).with_name('navigator.sqlite3'))
    for code in ('WOUND', 'ANGIO', 'NERVE'):
        print(code, [(x['parent'], x['variant']) for x in recommend(r,s,p,code)])
