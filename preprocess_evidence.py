"""Tissue Regeneration Peptide Navigator: conservative evidence preprocessing.

Reads existing extraction_candidates and papers; writes ONLY processed_evidence.
No source Excel or existing evidence is modified. Does not verify efficacy.
Run: python preprocess_evidence.py --db navigator.sqlite3
     python preprocess_evidence.py --db navigator.sqlite3 --report
     python preprocess_evidence.py --self-test
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

DB_DEFAULT = Path(__file__).resolve().parent / 'navigator.sqlite3'
PEPTIDES = ('GHK', 'RADA16-I', 'QK', 'IKVAV', 'KLD-12', 'RGD', 'YIGSR')
MISSING = 'Not reported'

# Recognize only explicit, tightly constrained comparisons. All are candidates,
# NEVER evidence considered validated for the recommendation score.
NUM = r'(\d+(?:\.\d+)?)'
PCT = rf'{NUM}\s*%'
TREAT = r'(?:treated|treatment|experimental|peptide[- ]treated|modified[- ]peptide)\s+(?:group|cells|samples)?'
CTRL = r'(?:untreated|vehicle|control|placebo)\s+(?:group|cells|samples)?'
METRIC = re.compile(r'wound closure|scratch closure|cell migration|cell viability|tube formation|neurite (?:length|outgrowth)|axon (?:length|growth)|proliferation|collagen (?:deposition|production)|vascularization', re.I)
NEGATION = re.compile(r'\b(?:not|no|without|neither|failed to|did not|was not|were not)\b.{0,40}\b(?:increase|enhance|improve|promote|decrease|reduce|significant)', re.I)
DIRECTION = re.compile(r'\b(increased?|enhanced?|improved?|promoted?|decreased?|reduced?|inhibited?|higher|lower)\b', re.I)
PAIRED = [
    re.compile(rf'\b{TREAT}\b[^.;]{{0,90}}?{PCT}[^.;]{{0,130}}?\b(?:compared (?:with|to)|versus|vs\.?)\b[^.;]{{0,90}}?\b{CTRL}\b[^.;]{{0,75}}?{PCT}', re.I),
    re.compile(rf'\b{PCT}[^.;]{{0,55}}?\b(?:in|for|of)\s+(?:the\s+)?{TREAT}\b[^.;]{{0,100}}?\b(?:compared (?:with|to)|versus|vs\.?)\b[^.;]{{0,70}}?{PCT}[^.;]{{0,45}}?\b(?:in|for|of)\s+(?:the\s+)?{CTRL}\b', re.I),
]


def now():
    return datetime.now(timezone.utc).isoformat()


def clean(text):
    return re.sub(r'\s+', ' ', str(text or '')).strip()


def doi_normalized(text):
    return re.sub(r'^https?://(?:dx\.)?doi\.org/', '', clean(text), flags=re.I).lower()


def explicitly_mentioned_peptide(text, peptide):
    aliases = {
        'GHK': r'(?<![A-Za-z])(?:GHK(?:-Cu|-C16)?|Cu-GHK(?:-C16)?)(?![A-Za-z])',
        'RADA16-I': r'\bRADA16(?:-I)?\b',
        'QK': r'\bQK\s+peptide\b',
        'IKVAV': r'\bIKVAV\b',
        'KLD-12': r'\bKLD[- ]?12\b',
        'RGD': r'\bRGD\s+(?:peptide|motif|sequence)\b',
        'YIGSR': r'\bYIGSR\b',
    }
    return bool(re.search(aliases[peptide], text, re.I))


def variant_mention(text, peptide):
    if peptide == 'GHK':
        for label, pattern in (
            ('Cu-GHK-C16', r'\b(?:Cu-GHK-C16|GHK-C16-Cu)\b'),
            ('GHK-C16', r'\bGHK-C16\b'),
            ('GHK-Cu', r'\b(?:GHK-Cu|Cu-GHK)\b'),
        ):
            if re.search(pattern, text, re.I):
                return label
    return peptide if explicitly_mentioned_peptide(text, peptide) else MISSING


def classify(text, peptide):
    """Conservative extraction. Percent calculations only for explicit pairs.

    Even valid mathematical parsing remains AUTO_EXTRACTED, not verified efficacy.
    """
    text = clean(text)
    metric = METRIC.search(text)
    named = explicitly_mentioned_peptide(text, peptide)
    variant = variant_mention(text, peptide)
    base = {
        'matched_peptide': peptide,
        'matched_variant': variant,
        'metric': metric.group(0) if metric else MISSING,
        'treatment_value': None,
        'control_value': None,
        'unit': MISSING,
        'absolute_difference_pp': None,
        'relative_change_percent': None,
        'reported_relative_change_percent': None,
        'comparison': MISSING,
        'source_sentence': text,
        'qualitative_statement': f'원문에 “{text}”라고 명시됨',
        'extraction_status': 'REVIEW_REQUIRED',
        'review_reason': '',
        'verified_for_ranking': False,
    }
    if not named:
        base['review_reason'] = '문장 안에 대상 펩타이드 명칭이 직접 확인되지 않음 (논문 차원의 연결만 존재)'
        return base
    if NEGATION.search(text):
        base['review_reason'] = '부정/비유의 문맥: 양성 효과로 해석하거나 수치 비교하지 않음'
        return base
    if not metric:
        base['review_reason'] = '평가 지표를 특정할 수 없음'
        return base
    if re.search(r'\b(?:respectively|multiple groups|different concentrations)\b', text, re.I):
        base['review_reason'] = '복수 실험군/조건 가능성: 비교 수치 연결 불명확'
        return base
    for pattern in PAIRED:
        match = pattern.search(text)
        if not match:
            continue
        # Two explicitly tied percentages must be the ONLY percentages in the sentence.
        if len(re.findall(PCT, text)) != 2:
            base['review_reason'] = '숫자가 둘 이상 섞여 있어 실험군-대조군 확정 불가'
            return base
        treat, ctrl = float(match.group(1)), float(match.group(2))
        if not all(0 <= v <= 100 for v in (treat, ctrl)):
            base['review_reason'] = '백분율 범위 오류'
            return base
        base.update(treatment_value=treat, control_value=ctrl, unit='%',
                    absolute_difference_pp=round(treat - ctrl, 4),
                    relative_change_percent=round((treat - ctrl) / ctrl * 100, 4) if ctrl else None,
                    comparison='명시된 treatment vs control (세부 조건 원문 확인 필요)',
                    extraction_status='AUTO_EXTRACTED_PAIRED_VALUES',
                    review_reason='백분율 쌍을 자동 추출함. 실험 조건, 측정 시점, 문장 분할 및 군 연결을 원저로 확인해야 함')
        if ctrl == 0:
            base['review_reason'] += '; control=0이라 상대 변화율 계산 불가'
        return base
    # "increased by 45% compared with control" is a reported relative change,
    # NOT a reconstructed pair of raw treatment and control values.
    explicit = re.search(r'\b(increased?|decreased?|reduced?)\b[^.;]{0,70}?\bby\s+(\d+(?:\.\d+)?)\s*%\s+(?:compared\s+(?:with|to)|relative\s+to)\s+(?:the\s+)?(?:untreated\s+|vehicle\s+)?control\b', text, re.I)
    if explicit and len(re.findall(PCT, text)) == 1:
        val = float(explicit.group(2))
        if val > 10000:
            base['review_reason'] = '상대 변화율 범위 검토 필요'
            return base
        base['reported_relative_change_percent'] = -val if re.match(r'decreas|reduc', explicit.group(1), re.I) else val
        base['comparison'] = '문장에 명시된 control 대비 변화 (세부 군 확인 필요)'
        base['extraction_status'] = 'AUTO_EXTRACTED_REPORTED_CHANGE'
        base['review_reason'] = '원문 표현의 변화율을 기록한 것으로 직접 계산한 값이 아님. 원저 확인 필요'
        return base
    base['extraction_status'] = 'QUALITATIVE_QUOTE'
    base['review_reason'] = '명확히 연결된 treatment/control 정량값이 없어 원문 인용으로 보존'
    return base


def initialize(con):
    needed = {'extraction_candidates', 'papers'}
    existing = {r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if not needed.issubset(existing):
        raise ValueError('기존 navigator.sqlite3 구조가 아님: ' + ', '.join(sorted(needed-existing)))
    con.execute('''CREATE TABLE IF NOT EXISTS processed_evidence (
      candidate_id TEXT PRIMARY KEY,
      paper_key TEXT NOT NULL,
      peptide TEXT NOT NULL,
      function_code TEXT,
      source_section TEXT,
      doi TEXT,
      paper_url TEXT,
      source_text_hash TEXT NOT NULL,
      processor_version TEXT NOT NULL,
      result_json TEXT NOT NULL,
      extraction_status TEXT NOT NULL,
      processed_at TEXT NOT NULL
    )''')
    con.execute('CREATE INDEX IF NOT EXISTS idx_processed_function ON processed_evidence(function_code, peptide)')
    con.commit()


def process(con, force=False):
    data = con.execute('''SELECT c.candidate_id,c.paper_key,c.peptide,c.source_section,c.source_text,
                                 c.screening_json,p.doi,p.paper_url
                          FROM extraction_candidates c JOIN papers p ON c.paper_key=p.paper_key
                          ORDER BY c.candidate_id''').fetchall()
    changes = 0
    for cid, key, pep, section, text, screening_json, paper_doi, paper_url in data:
        if pep not in PEPTIDES:
            continue
        h = hashlib.sha256(text.encode('utf-8')).hexdigest()
        old = con.execute('SELECT source_text_hash,processor_version FROM processed_evidence WHERE candidate_id=?', (cid,)).fetchone()
        if not force and old == (h, '1.0'):
            continue
        try:
            screening = json.loads(screening_json or '{}')
        except json.JSONDecodeError:
            screening = {}
        result = classify(text, pep)
        function = clean(screening.get('Target_Function_Code'))
        con.execute('''INSERT INTO processed_evidence
          (candidate_id,paper_key,peptide,function_code,source_section,doi,paper_url,
           source_text_hash,processor_version,result_json,extraction_status,processed_at)
          VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
          ON CONFLICT(candidate_id) DO UPDATE SET
          paper_key=excluded.paper_key,peptide=excluded.peptide,function_code=excluded.function_code,
          source_section=excluded.source_section,doi=excluded.doi,paper_url=excluded.paper_url,
          source_text_hash=excluded.source_text_hash,processor_version=excluded.processor_version,
          result_json=excluded.result_json,extraction_status=excluded.extraction_status,
          processed_at=excluded.processed_at''',
          (cid,key,pep,function,section,doi_normalized(paper_doi),paper_url,h,'1.0',
           json.dumps(result,ensure_ascii=False),result['extraction_status'],now()))
        changes += 1
    con.commit()
    return changes


def report(con):
    total = con.execute('SELECT COUNT(*) FROM processed_evidence').fetchone()[0]
    print(f'processed_evidence {total}')
    for name, n in con.execute('SELECT extraction_status, COUNT(*) FROM processed_evidence GROUP BY extraction_status ORDER BY extraction_status'):
        print(name, n)
    print('NOTE: 자동 전처리 결과는 검증 완료 근거가 아니며 추천 점수를 높이지 않습니다.')


def self_test():
    examples = [
        ('GHK increased wound closure in the treatment group to 60% compared with the control group at 40%.', 'GHK', 'AUTO_EXTRACTED_PAIRED_VALUES', 50.0),
        ('GHK enhanced cell migration compared with the untreated control group.', 'GHK', 'QUALITATIVE_QUOTE', None),
        ('GHK increased wound closure by 30% compared with control.', 'GHK', 'AUTO_EXTRACTED_REPORTED_CHANGE', None),
        ('GHK did not significantly improve wound closure compared with control.', 'GHK', 'REVIEW_REQUIRED', None),
        ('The treatment group showed 60% wound closure versus 40% in control.', 'GHK', 'REVIEW_REQUIRED', None),
    ]
    for sentence, peptide, status, relative in examples:
        got = classify(sentence,peptide)
        assert got['extraction_status'] == status, (status, got)
        assert got['relative_change_percent'] == relative, got
        assert not got['verified_for_ranking']
    print('SELF-TEST PASS (5 synthetic sentences; no original paper verification)')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--db',type=Path,default=DB_DEFAULT)
    parser.add_argument('--force',action='store_true',help='전처리 결과 재생성 (원본 후보 보존)')
    parser.add_argument('--report',action='store_true',help='상태 조회만')
    parser.add_argument('--self-test',action='store_true')
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if not args.db.exists():
        parser.error(f'SQLite 파일이 없습니다: {args.db}')
    with sqlite3.connect(args.db) as con:
        initialize(con)
        if not args.report:
            print('NEW_OR_UPDATED',process(con,args.force))
        report(con)

if __name__ == '__main__':
    main()
