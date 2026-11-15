"""Evidence preprocessing v3 for Tissue Regeneration Peptide Navigator.

Reads existing extraction_candidates and optional builder v2 context. Never alters
source papers or extraction candidates. Quantitative matches are UNVERIFIED until
human review; figure-only values are never guessed.

Commands:
  python preprocess_evidence.py --self-test
  python preprocess_evidence.py --db navigator.sqlite3
  python preprocess_evidence.py --db navigator.sqlite3 --report
  python preprocess_evidence.py --db navigator.sqlite3 --force
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_DB = Path(__file__).resolve().parent / 'navigator.sqlite3'
VERSION = '3.0'
PEPTIDES = ('GHK', 'RADA16-I', 'QK', 'IKVAV', 'KLD-12', 'RGD', 'YIGSR')
NA = 'Not reported'

ALIASES = {
    'GHK': r'(?<![A-Za-z])(?:GHK(?:[- ](?:Cu|C16))?|Cu[- ]GHK(?:[- ]C16)?|glycyl[- ]histidyl[- ]lysine)(?![A-Za-z])',
    'RADA16-I': r'(?<![A-Za-z])RADA16(?:[- ](?:I|RGD|IKVAV))?(?![A-Za-z])',
    'QK': r'\b(?:QK\s+peptide|QK[- ](?:GelMA|DA)|VEGF mimetic peptide)\b',
    'IKVAV': r'\bIKVAV\b',
    'KLD-12': r'\bKLD[- ]?12\b',
    'RGD': r'\b(?:RGD\s+(?:peptide|motif|sequence)|RADA16[- ]RGD)\b',
    'YIGSR': r'\bYIGSR\b',
}
METRICS = {
    'Wound closure': r'wound (?:closure|area|contraction|healing)|re[- ]?epithelialization|scratch (?:assay|closure)',
    'Cell migration': r'(?:cell|keratinocyte|fibroblast) migration|migrat(?:ion|ed|ing)',
    'Cell proliferation': r'proliferat(?:ion|ed|ing)|cell growth',
    'Cell viability': r'cell viability|cell survival|cytotoxicity',
    'Angiogenesis': r'angiogen(?:esis|ic)|vasculari[sz]ation|tube formation|CD31[- ]positive|blood vessel',
    'Axon/neurite': r'axon(?:al)? (?:growth|regenerat\w*|extension|outgrowth|length)|neurite (?:outgrowth|length|extension)|nerve regenerat\w*',
    'Schwann cells': r'Schwann cell(?:s)? (?:infiltrat\w*|migration|proliferat\w*)',
    'Collagen/ECM': r'collagen (?:synthesis|production|deposition)|extracellular matrix (?:deposition|remodeling)',
}
NEGATIVE = re.compile(r'\b(?:did not|does not|do not|was not|were not|no significant|not significant|non[- ]significant|failed to|no difference)\b', re.I)
DIRECTION = re.compile(r'\b(?:increase\w*|decrease\w*|enhanc\w*|improv\w*|promot\w*|reduc\w*|inhibit\w*|higher|lower|upregulat\w*|downregulat\w*|attenuat\w*|accelerat\w*|greater|smaller|more|less|fold)\b', re.I)
CITATION = re.compile(r'(?:\[\s*\d+(?:\s*[,–-]\s*\d+)*\s*\]|\(\s*[A-Z][A-Za-z-]+\s+et al\.?[, ]+\d{4}\s*\)|\b[A-Z][a-z]+ et al\.?)')
GROUP_T = r'(?:treated|treatment|experimental|peptide[- ]treated|intervention)'
GROUP_C = r'(?:untreated|control|vehicle|blank)'
NUM = r'(?P<{k}>\d+(?:\.\d+)?)'
# Explicitly-labelled same-unit percentages in one sentence; avoid silent ordering assumptions.
PAIR_PATTERNS = [
    re.compile(rf'{NUM.format(k="t")}\s*%\s*(?:in|for|of|among)\s+(?:the\s+)?{GROUP_T}\s*(?:group|cells)?[^.;]{{0,100}}?(?:compared\s+(?:with|to)|versus|vs\.?|and)[^.;]{{0,60}}?{NUM.format(k="c")}\s*%\s*(?:in|for|of|among)?\s*(?:the\s+)?{GROUP_C}\b', re.I),
    re.compile(rf'{GROUP_T}\s*(?:group|cells)?[^.;]{{0,65}}?{NUM.format(k="t")}\s*%[^.;]{{0,100}}?(?:compared\s+(?:with|to)|versus|vs\.?|and)\s*(?:the\s+)?{GROUP_C}\s*(?:group|cells)?[^.;]{{0,65}}?{NUM.format(k="c")}\s*%', re.I),
    re.compile(rf'{GROUP_C}\s*(?:group|cells)?[^.;]{{0,65}}?{NUM.format(k="c")}\s*%[^.;]{{0,100}}?(?:compared\s+(?:with|to)|versus|vs\.?|and)\s*(?:the\s+)?{GROUP_T}\s*(?:group|cells)?[^.;]{{0,65}}?{NUM.format(k="t")}\s*%', re.I),
]
REPORTED = re.compile(r'\b(?P<dir>increas\w*|decreas\w*|reduc\w*|improv\w*|enhanc\w*)\b[^.;]{0,70}?\bby\s+(?P<n>\d+(?:\.\d+)?)\s*%\s+(?:compared\s+(?:with|to)|relative\s+to|versus|vs\.?)\s+(?:the\s+)?(?:untreated\s+|vehicle\s+)?control\b', re.I)
FOLD = re.compile(r'\b(?P<n>\d+(?:\.\d+)?)\s*[- ]?fold\s+(?P<dir>increase|decrease|higher|lower)\b[^.;]{0,90}?(?:compared\s+(?:with|to)|relative\s+to|versus|vs\.?)\s+(?:the\s+)?(?:untreated\s+|vehicle\s+)?control\b', re.I)
STATS = re.compile(r'(?<!\w)p\s*(?:<|≤|=|>|≥)\s*0?\.\d+', re.I)
MEAN_SD = re.compile(r'(?P<mean>\d+(?:\.\d+)?)\s*(?:±|\+/-)\s*(?P<spread>\d+(?:\.\d+)?)\s*(?P<unit>%|mm|µm|μm|um|cm|cells?|a\.?u\.?)?', re.I)
FIGURE = re.compile(r'\b(?:fig(?:ure)?\.?\s*\d+[A-Za-z]?|table\s*\d+)\b', re.I)
MULTI = re.compile(r'\b(?:respectively|different concentrations|various doses|multiple time points)\b', re.I)
REVIEW_TITLE = re.compile(r'\breview\b|meta[- ]analysis|bibliometric|systematic review|scoping review', re.I)


def clean(x):
    return re.sub(r'\s+', ' ', str(x or '')).strip()


def norm_doi(x):
    return re.sub(r'^https?://(?:dx\.)?doi\.org/', '', clean(x), flags=re.I).lower()


def variant(text, peptide):
    pats = {
        'GHK': [('Cu-GHK-C16', r'\b(?:Cu-GHK-C16|GHK-C16-Cu)\b'),
                ('GHK-C16', r'\bGHK-C16\b'), ('GHK-Cu', r'\b(?:GHK-Cu|Cu-GHK)\b')],
        'RADA16-I': [('RADA16-RGD', r'\bRADA16-RGD\b'),
                     ('RADA16-IKVAV', r'\bRADA16-IKVAV\b')],
    }
    for label, pattern in pats.get(peptide, []):
        if re.search(pattern, text, re.I):
            return label
    return peptide


def source_risk(title, study_type, has_citation, is_primary, section):
    if study_type in ('review', 'review_possible') or REVIEW_TITLE.search(title):
        return '리뷰 논문: 해당 논문의 독립 실험 결과로 사용 불가'
    if has_citation:
        return '인용 표시가 포함됨: 타 연구 결과의 재인용 여부 확인 필요'
    if study_type == 'unknown':
        return '원저 여부 미확인: 자료는 추출하되 직접 근거로 확정 불가'
    if is_primary is False:
        return '직접 실험 근거로 표시되지 않은 후보'
    if section and re.search(r'\bdiscussion\b', section, re.I) and not re.search(r'\bresults?\b', section, re.I):
        return 'Discussion 문장: 원문 Results 재확인 필요'
    return ''


def extract(sentence, peptide, title='', section='', context=None, study_type=None):
    """Return structured but UNVERIFIED data. Context cannot silently supply group values."""
    sentence = clean(sentence)
    context = context or {}
    paragraph = clean(context.get('paragraph'))
    previous = clean(context.get('prior_paragraph'))
    # An exact peptide mention in the same source sentence is required for attribution.
    has_peptide = bool(re.search(ALIASES[peptide], sentence, re.I))
    metrics = [k for k, p in METRICS.items() if re.search(p, sentence, re.I)]
    contextual_metrics = [k for k, p in METRICS.items() if re.search(p, paragraph, re.I)]
    cited = bool(CITATION.search(sentence)) or bool(context.get('citation_warning'))
    primary = context.get('is_primary_candidate')
    status = 'REVIEW_REQUIRED'
    reason = ''
    out = {
        'matched_peptide': peptide,
        'matched_variant': variant(sentence, peptide) if has_peptide else NA,
        'metric': '; '.join(metrics) or NA,
        'contextual_metrics': '; '.join(contextual_metrics) or NA,
        'treatment_value': None, 'control_value': None,
        'unit': NA, 'absolute_difference_pp': None,
        'relative_change_percent': None, 'reported_relative_change_percent': None,
        'reported_fold_change': None, 'mean_spread_values': [],
        'p_values': STATS.findall(sentence),
        'figure_or_table_mentions': FIGURE.findall(sentence),
        'comparison': NA, 'source_sentence': sentence,
        'qualitative_statement': f'원문에 “{sentence}”라고 명시됨',
        'source_section': section,
        'source_paragraph': paragraph, 'prior_paragraph': previous,
        'context_link_status': 'sentence_explicit' if has_peptide else 'unresolved',
        'has_citation_marker': cited,
        'study_type': study_type or 'unknown',
        'verified_for_ranking': False,
        'extraction_status': status, 'review_reason': '',
        'review_article_possible': bool(REVIEW_TITLE.search(title)),
        'quantitative_source': NA,
    }
    for m in MEAN_SD.finditer(sentence):
        out['mean_spread_values'].append({
            'mean': float(m.group('mean')), 'spread': float(m.group('spread')),
            'unit': m.group('unit') or NA,
            'spread_type': 'SD_or_SE_or_other_unresolved',
        })

    risk = source_risk(title, study_type, cited, primary, section)
    if risk:
        out['review_reason'] = risk
        return out
    if not has_peptide:
        out['review_reason'] = '해당 문장에서 대상 펩타이드 확인 불가; 문단 문맥만으로 귀속 불가'
        return out
    if not metrics:
        out['review_reason'] = '결과 지표가 같은 문장에 명확하지 않음; 문단 지표는 참고만 저장'
        return out
    if NEGATIVE.search(sentence):
        out['review_reason'] = '부정 또는 비유의 결과: 원문 보존, 우수 효능으로 자동 판단 불가'
        out['effect_direction'] = 'no_effect_or_negative_reported'
        return out
    if not DIRECTION.search(sentence):
        out['review_reason'] = '효과 방향 표현 미확인'
        return out
    if MULTI.search(sentence):
        out['review_reason'] = '복수 군·시간·농도 표기: 수치쌍의 조건 연결 확인 필요'
        return out

    percentages = re.findall(r'\d+(?:\.\d+)?\s*%', sentence)
    for pattern in PAIR_PATTERNS:
        m = pattern.search(sentence)
        if m and len(percentages) == 2:
            t, c = float(m.group('t')), float(m.group('c'))
            if 0 <= t <= 100 and 0 <= c <= 100:
                out.update(
                    treatment_value=t, control_value=c, unit='%',
                    absolute_difference_pp=round(t-c, 4),
                    relative_change_percent=round((t-c)/c*100, 4) if c != 0 else None,
                    comparison='explicit treatment vs control; time/dose validation pending',
                    quantitative_source='same_sentence_group_pair',
                    extraction_status='AUTO_EXTRACTED_PAIRED_VALUES',
                    review_reason='명시된 동일 단위 수치쌍; 실험 조건 및 원문 검증 필요' + ('; control=0' if c == 0 else ''),
                )
                return out

    m = REPORTED.search(sentence)
    if m and len(percentages) == 1:
        amount = float(m.group('n'))
        sign = -1 if re.match(r'decreas|reduc', m.group('dir'), re.I) else 1
        out.update(
            reported_relative_change_percent=round(sign*amount, 4),
            comparison='reported relative change against control',
            quantitative_source='author_reported_percent',
            extraction_status='AUTO_EXTRACTED_REPORTED_CHANGE',
            review_reason='원문 보고 상대 변화율; 처리 조건과 방향 검증 필요',
        )
        return out

    m = FOLD.search(sentence)
    if m:
        amount = float(m.group('n'))
        direction = m.group('dir').lower()
        # '1.8-fold decrease' does not unambiguously encode treatment/control
        # ratio across conventions. Save literal only; do not invent percentage.
        out.update(
            reported_fold_change=amount,
            comparison=f'author-reported {direction} relative to control',
            quantitative_source='author_reported_fold',
            extraction_status='AUTO_EXTRACTED_FOLD',
            review_reason='원문 fold 값 저장; fold decrease 등의 해석 및 비교 조건 검증 필요',
        )
        return out

    if out['mean_spread_values']:
        out.update(
            quantitative_source='unassigned_mean_spread',
            extraction_status='QUANTITATIVE_UNPAIRED',
            review_reason='평균±분산값 발견; 실험군·대조군 및 SD/SEM 구분 필요',
        )
        return out

    if out['p_values']:
        out.update(
            quantitative_source='statistical_only',
            extraction_status='QUANTITATIVE_UNPAIRED',
            review_reason='통계값 발견; 비교군과 효과 크기는 미확인',
        )
        return out

    if out['figure_or_table_mentions']:
        out.update(
            extraction_status='FIGURE_TABLE_REVIEW',
            review_reason='Figure/Table에 정량 데이터가 존재할 수 있음; 그래프·표 수치는 자동 추정하지 않음',
        )
        return out

    out.update(
        extraction_status='QUALITATIVE_QUOTE',
        review_reason='효과 관련 원문 표현 발견; 정량 비교는 확인되지 않음',
    )
    return out


def initialize(db):
    tables = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if not {'papers', 'extraction_candidates'}.issubset(tables):
        raise ValueError('필수 테이블 papers/extraction_candidates가 없습니다')
    db.execute('''CREATE TABLE IF NOT EXISTS processed_evidence (
      candidate_id TEXT PRIMARY KEY, paper_key TEXT NOT NULL, peptide TEXT NOT NULL,
      function_code TEXT, source_section TEXT, doi TEXT, paper_url TEXT,
      source_text_hash TEXT NOT NULL, processor_version TEXT NOT NULL,
      result_json TEXT NOT NULL, extraction_status TEXT NOT NULL, processed_at TEXT NOT NULL)''')
    db.execute('CREATE INDEX IF NOT EXISTS idx_processed_function ON processed_evidence(function_code,peptide)')
    db.execute('''CREATE TABLE IF NOT EXISTS evidence_audit_v3 (
      candidate_id TEXT PRIMARY KEY, paper_key TEXT NOT NULL, peptide TEXT NOT NULL,
      study_type TEXT NOT NULL, has_context INTEGER NOT NULL,
      quantitative_source TEXT NOT NULL, verification_status TEXT NOT NULL,
      audit_json TEXT NOT NULL, updated_at TEXT NOT NULL)''')
    db.commit()
    return tables


def process(db, force=False):
    existing = {row[0] for row in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    has_context = 'candidate_context_v2' in existing
    has_screen = 'paper_screening_v2' in existing
    context_join = '''LEFT JOIN candidate_context_v2 x ON c.candidate_id=x.candidate_id''' if has_context else ''
    screen_join = '''LEFT JOIN paper_screening_v2 s ON c.paper_key=s.paper_key''' if has_screen else ''
    context_select = 'x.paragraph,x.prior_paragraph,x.citation_warning,x.is_primary_candidate' if has_context else 'NULL,NULL,NULL,NULL'
    screen_select = 's.study_type' if has_screen else 'NULL'
    rows = db.execute(f'''
      SELECT c.candidate_id,c.paper_key,c.peptide,c.source_section,c.source_text,
        c.screening_json,p.doi,p.paper_url,p.title,
        {context_select},{screen_select}
      FROM extraction_candidates c
      JOIN papers p ON c.paper_key=p.paper_key
      {context_join} {screen_join}
      ORDER BY c.candidate_id''').fetchall()
    changed = 0
    for cid, key, peptide, sec, sentence, raw, doi, url, title, para, previous, citation_flag, primary, study_type in rows:
        if peptide not in PEPTIDES:
            continue
        context = {
            'paragraph': para or '', 'prior_paragraph': previous or '',
            'citation_warning': citation_flag,
            'is_primary_candidate': bool(primary) if primary is not None else None,
        }
        fingerprint = hashlib.sha256(json.dumps(
            [sentence,title,sec,context,study_type], ensure_ascii=False,
            sort_keys=True).encode()).hexdigest()
        old = db.execute('SELECT source_text_hash,processor_version FROM processed_evidence WHERE candidate_id=?',(cid,)).fetchone()
        if not force and old == (fingerprint, VERSION):
            continue
        try:
            payload = json.loads(raw or '{}')
        except (TypeError, ValueError):
            payload = {}
        result = extract(sentence, peptide, title, sec, context, study_type)
        timestamp = datetime.now(timezone.utc).isoformat()
        db.execute('''INSERT INTO processed_evidence VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
          ON CONFLICT(candidate_id) DO UPDATE SET
            paper_key=excluded.paper_key,peptide=excluded.peptide,
            function_code=excluded.function_code,source_section=excluded.source_section,
            doi=excluded.doi,paper_url=excluded.paper_url,
            source_text_hash=excluded.source_text_hash,
            processor_version=excluded.processor_version,result_json=excluded.result_json,
            extraction_status=excluded.extraction_status,processed_at=excluded.processed_at''',
            (cid,key,peptide,clean(payload.get('Target_Function_Code')),
             sec,norm_doi(doi),url,fingerprint,VERSION,json.dumps(result,ensure_ascii=False),
             result['extraction_status'],timestamp))
        audit = {
            'candidate_id': cid,
            'reason': result['review_reason'],
            'source_section': sec,
            'doi': norm_doi(doi),
            'paper_url': url,
            'manual_checks': [
                'source publication type', 'correct intervention and comparator',
                'same assay, dose and time', 'quantitative value matches original source',
                'statistical support and conflicting outcomes',
            ],
        }
        db.execute('''INSERT INTO evidence_audit_v3 VALUES (?,?,?,?,?,?,?,?,?)
          ON CONFLICT(candidate_id) DO UPDATE SET
            paper_key=excluded.paper_key,peptide=excluded.peptide,
            study_type=excluded.study_type,has_context=excluded.has_context,
            quantitative_source=excluded.quantitative_source,
            verification_status=excluded.verification_status,
            audit_json=excluded.audit_json,updated_at=excluded.updated_at''',
            (cid,key,peptide,study_type or 'unknown',int(bool(para)),
             result['quantitative_source'],'UNVERIFIED',json.dumps(audit,ensure_ascii=False),timestamp))
        changed += 1
    db.commit()
    return changed


def report(db):
    n = db.execute('SELECT COUNT(*) FROM processed_evidence').fetchone()[0]
    print('processed_evidence',n)
    for name, count in db.execute('SELECT extraction_status,COUNT(*) FROM processed_evidence GROUP BY extraction_status ORDER BY extraction_status'):
        print(name,count)
    print('=== 검토 사유 TOP 8 ===')
    reasons=Counter()
    for (raw,) in db.execute('SELECT result_json FROM processed_evidence'):
        try:
            item=json.loads(raw)
            reasons[item.get('review_reason') or '미기재']+=1
        except (ValueError,TypeError):
            pass
    for reason,count in reasons.most_common(8):
        print(count,reason)
    for table in ('candidate_context_v2','paper_screening_v2','evidence_audit_v3'):
        if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?",(table,)).fetchone():
            print(table,db.execute(f'SELECT COUNT(*) FROM {table}').fetchone()[0])
    print('주의: 추출 결과는 모두 UNVERIFIED이며 검증 없이 후보 추천 점수에 사용하지 않습니다.')
    print('현재 builder v2는 Table의 셀 값이나 Figure 그래프 픽셀을 수집하지 않습니다.')


def self_test():
    samples = [
        ('GHK increased wound closure to 60% in treated cells compared with 40% in control cells.', 'GHK', 'AUTO_EXTRACTED_PAIRED_VALUES', 50.0),
        ('GHK increased wound closure by 30% compared with control.', 'GHK', 'AUTO_EXTRACTED_REPORTED_CHANGE', None),
        ('GHK enhanced cell migration 1.8-fold higher compared with control.', 'GHK', 'AUTO_EXTRACTED_FOLD', None),
        ('GHK increased cell migration (23.4 ± 3.1%) after treatment.', 'GHK', 'QUANTITATIVE_UNPAIRED', None),
        ('GHK improved wound closure (p < 0.05).', 'GHK', 'QUANTITATIVE_UNPAIRED', None),
        ('GHK promoted axonal regeneration (Figure 2B).', 'GHK', 'FIGURE_TABLE_REVIEW', None),
        ('GHK promoted cell migration.', 'GHK', 'QUALITATIVE_QUOTE', None),
        ('GHK did not improve wound closure.', 'GHK', 'REVIEW_REQUIRED', None),
        ('GHK not only improved biosafety but also enhanced collagen synthesis.', 'GHK', 'QUALITATIVE_QUOTE', None),
        ('RADA16-RGD promoted axonal regeneration [172].', 'RADA16-I', 'REVIEW_REQUIRED', None),
        ('Other materials promoted cell migration.', 'GHK', 'REVIEW_REQUIRED', None),
    ]
    for sentence,peptide,status,change in samples:
        output=extract(sentence,peptide,study_type='original_candidate')
        assert output['extraction_status']==status, (sentence,output)
        assert output['relative_change_percent']==change, (sentence,output)
        assert not output['verified_for_ranking']
    review=extract('GHK enhanced cell migration.', 'GHK', study_type='review')
    assert review['extraction_status']=='REVIEW_REQUIRED'
    ambiguous=extract('GHK increased wound closure in treated cells and control cells respectively.', 'GHK',study_type='original_candidate')
    assert ambiguous['extraction_status']=='REVIEW_REQUIRED'
    print(f'SELF-TEST PASS: {len(samples)+2} synthetic checks; no paper-level accuracy claimed')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--db',type=Path,default=DEFAULT_DB)
    ap.add_argument('--self-test',action='store_true')
    ap.add_argument('--report',action='store_true')
    ap.add_argument('--force',action='store_true')
    args=ap.parse_args()
    if args.self_test:
        self_test(); return
    if not args.db.is_file():
        ap.error('DB를 찾을 수 없습니다: '+str(args.db))
    with sqlite3.connect(args.db) as db:
        initialize(db)
        if not args.report:
            print('NEW_OR_UPDATED',process(db,args.force))
        report(db)


if __name__=='__main__':
    main()

