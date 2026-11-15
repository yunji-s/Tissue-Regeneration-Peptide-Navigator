"""Conservative literature preprocessing for Tissue Regeneration Peptide Navigator.

Writes processed_evidence only; preserves original source rows. Run:
 python preprocess_evidence.py --self-test
 python preprocess_evidence.py --db navigator.sqlite3
 python preprocess_evidence.py --db navigator.sqlite3 --report
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

ROOT = Path(__file__).resolve().parent
DEFAULT_DB = ROOT / 'navigator.sqlite3'
VERSION = '2.0'
PEPTIDES = ('GHK', 'RADA16-I', 'QK', 'IKVAV', 'KLD-12', 'RGD', 'YIGSR')
ALIASES = {
    'GHK': r'(?<![A-Za-z])(?:GHK(?:[- ](?:Cu|C16))?|Cu[- ]GHK(?:[- ]C16)?)(?![A-Za-z])',
    'RADA16-I': r'(?<![A-Za-z])RADA16(?:[- ]I|[- ]RGD|[- ]IKVAV)?(?![A-Za-z])',
    'QK': r'\bQK\s+peptide\b',
    'IKVAV': r'\bIKVAV\b',
    'KLD-12': r'\bKLD[- ]?12\b',
    'RGD': r'\bRGD(?:\s+(?:peptide|motif|sequence))?\b',
    'YIGSR': r'\bYIGSR\b',
}
METRICS = {
    'Wound closure': r'wound (?:closure|area|contraction|healing)|re[- ]?epithelialization|scratch (?:assay|closure)',
    'Cell migration': r'(?:cell|keratinocyte|fibroblast) migration|migrat(?:ion|ed|ing)',
    'Cell proliferation': r'proliferat(?:ion|ed|ing)|cell growth',
    'Cell viability': r'cell viability|cell survival|cytotoxicity',
    'Angiogenesis': r'angiogen(?:esis|ic)|vasculari[sz]ation|tube formation|CD31[- ]positive|blood vessel',
    'Axon/neurite': r'axon(?:al)? (?:growth|regenerat\w*|extension|outgrowth|length)|neurite (?:outgrowth|length|extension)|nerve regenerat\w*',
    'Schwann cells': r'Schwann cell (?:infiltrat\w*|migration|proliferat\w*)',
    'Collagen/ECM': r'collagen (?:synthesis|production|deposition)|extracellular matrix (?:deposition|remodeling)',
}
NEGATIVE = re.compile(r'\b(?:did not|does not|do not|was not|were not|no significant|not significant|non[- ]significant|failed to)\b', re.I)
EFFECT = re.compile(r'\b(?:increase\w*|decrease\w*|enhanc\w*|improv\w*|promot\w*|reduc\w*|inhibit\w*|higher|lower|upregulat\w*|downregulat\w*)\b',re.I)
CITATION = re.compile(r'(?:\[\s*\d+(?:\s*[,–-]\s*\d+)*\s*\]|\(\s*[A-Z][A-Za-z-]+\s+et al\.?[, ]+\d{4}\s*\)|\b[A-Z][a-z]+ et al\.?\s*\[?\d+\]?)')
PAIR = re.compile(r'(?P<t>\d+(?:\.\d+)?)\s*%\s*(?:in|for|of)\s+(?:the\s+)?(?:treated|treatment|experimental)\s+(?:group|cells)?[^.;]{0,110}?(?:compared\s+(?:with|to)|versus|vs\.?)[^.;]{0,60}?(?P<c>\d+(?:\.\d+)?)\s*%\s*(?:in|for|of)?\s*(?:the\s+)?(?:untreated|control|vehicle)\b',re.I)
PAIR2 = re.compile(r'(?:treated|treatment|experimental)\s+(?:group|cells)?[^.;]{0,80}?(?P<t>\d+(?:\.\d+)?)\s*%[^.;]{0,100}?(?:compared\s+(?:with|to)|versus|vs\.?)\s+(?:the\s+)?(?:untreated|control|vehicle)\s+(?:group|cells)?[^.;]{0,60}?(?:at|of|was|to)?\s*(?P<c>\d+(?:\.\d+)?)\s*%',re.I)
REPORTED = re.compile(r'\b(?P<dir>increas\w*|decreas\w*|reduc\w*)\b[^.;]{0,65}?\bby\s+(?P<n>\d+(?:\.\d+)?)\s*%\s+(?:compared\s+(?:with|to)|relative\s+to)\s+(?:the\s+)?(?:untreated\s+|vehicle\s+)?control\b',re.I)


def clean(s):
    return re.sub(r'\s+', ' ', str(s or '')).strip()


def norm_doi(s):
    return re.sub(r'^https?://(?:dx\.)?doi\.org/', '', clean(s), flags=re.I).lower()


def variant(text, peptide):
    if peptide == 'GHK':
        patterns = [('Cu-GHK-C16',r'\b(?:Cu-GHK-C16|GHK-C16-Cu)\b'),('GHK-C16',r'\bGHK-C16\b'),('GHK-Cu',r'\b(?:GHK-Cu|Cu-GHK)\b')]
    elif peptide == 'RADA16-I':
        patterns = [('RADA16-RGD',r'\bRADA16-RGD\b'),('RADA16-IKVAV',r'\bRADA16-IKVAV\b')]
    else:
        patterns = []
    for name, pat in patterns:
        if re.search(pat,text,re.I):
            return name
    return peptide


def classify(sentence, peptide, section='unknown', title=''):
    sentence = clean(sentence)
    names = re.findall(ALIASES[peptide],sentence,re.I) if peptide in ALIASES else []
    metrics = [name for name,pat in METRICS.items() if re.search(pat,sentence,re.I)]
    review_hint = bool(re.search(r'\breview\b|meta.analysis|systematic review|bibliometric',title,re.I))
    cited = bool(CITATION.search(sentence))
    result = dict(matched_peptide=peptide, matched_variant=variant(sentence,peptide) if names else 'Not reported',
                  metric='; '.join(metrics) or 'Not reported', treatment_value=None, control_value=None,
                  unit='Not reported', absolute_difference_pp=None, relative_change_percent=None,
                  reported_relative_change_percent=None, comparison='Not reported',
                  source_sentence=sentence, qualitative_statement=f'원문에 “{sentence}”라고 명시됨',
                  extraction_status='REVIEW_REQUIRED', review_reason='', verified_for_ranking=False,
                  source_section=section, review_article_possible=review_hint, has_citation_marker=cited)
    if review_hint or cited:
        result['review_reason']='리뷰 논문 또는 다른 연구 인용 가능성: 현재 논문의 직접 실험 근거로 사용 불가'
        return result
    if not names:
        result['review_reason']='해당 문장에 대상 펩타이드 명칭이 없음: 논문 연관성만으로 효능 연결 불가'
        return result
    if NEGATIVE.search(sentence):
        result['review_reason']='비유의/부정 결과: 원문 보존, 긍정적 효능 근거로 사용 불가'
        return result
    if not metrics:
        result['review_reason']='평가 지표 불명확: 추가 문맥 확인 필요'
        return result
    if not EFFECT.search(sentence):
        result['review_reason']='효과 방향 표현 미확인'
        return result
    if re.search(r'\brespectively\b|\bdifferent concentrations\b',sentence,re.I):
        result['review_reason']='복수 조건의 군별 수치 연결 불명확'
        return result
    values = re.findall(r'\d+(?:\.\d+)?\s*%',sentence)
    for pattern in (PAIR,PAIR2):
        m=pattern.search(sentence)
        if m and len(values)==2:
            treat,ctrl=float(m.group('t')),float(m.group('c'))
            if not (0<=treat<=100 and 0<=ctrl<=100):
                break
            result.update(treatment_value=treat,control_value=ctrl,unit='%',
                          absolute_difference_pp=round(treat-ctrl,3),
                          relative_change_percent=round((treat-ctrl)*100/ctrl,3) if ctrl else None,
                          comparison='명시된 treatment vs control (추가 검증 필요)',
                          extraction_status='AUTO_EXTRACTED_PAIRED_VALUES',
                          review_reason='자동 추출 수치쌍: 동일 실험·조건·측정 시점 여부 미검증'+('; 0으로 나눌 수 없음' if ctrl==0 else ''))
            return result
    m=REPORTED.search(sentence)
    if m and len(values)==1:
        n=float(m.group('n'))
        result['reported_relative_change_percent']=-n if re.match(r'decreas|reduc',m.group('dir'),re.I) else n
        result['comparison']='원문에 기재된 control 대비 변화'
        result['extraction_status']='AUTO_EXTRACTED_REPORTED_CHANGE'
        result['review_reason']='논문에 기재된 변화율로, 원자료 재계산값이 아님'
        return result
    result['extraction_status']='QUALITATIVE_QUOTE'
    result['review_reason']='연관 문장이지만 유효한 군별 정량 비교는 확인되지 않음'
    return result


def initialize(db):
    have={r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    if not {'extraction_candidates','papers'}.issubset(have):
        raise ValueError('기존 DB의 papers/extraction_candidates 테이블이 없습니다')
    db.execute('''CREATE TABLE IF NOT EXISTS processed_evidence (
        candidate_id TEXT PRIMARY KEY, paper_key TEXT NOT NULL, peptide TEXT NOT NULL,
        function_code TEXT, source_section TEXT, doi TEXT, paper_url TEXT,
        source_text_hash TEXT NOT NULL, processor_version TEXT NOT NULL,
        result_json TEXT NOT NULL, extraction_status TEXT NOT NULL, processed_at TEXT NOT NULL)''')
    db.execute('CREATE INDEX IF NOT EXISTS idx_processed_function ON processed_evidence(function_code,peptide)')
    db.commit()


def process(db, force=False):
    rows=db.execute('''SELECT c.candidate_id,c.paper_key,c.peptide,c.source_section,c.source_text,
        c.screening_json,p.doi,p.paper_url,p.title FROM extraction_candidates c
        JOIN papers p ON c.paper_key=p.paper_key ORDER BY c.candidate_id''').fetchall()
    count=0
    for cid,key,pep,sec,sentence,payload,doi,url,title in rows:
        if pep not in PEPTIDES:
            continue
        digest=hashlib.sha256((sentence+'\n'+clean(title)).encode()).hexdigest()
        old=db.execute('SELECT source_text_hash,processor_version FROM processed_evidence WHERE candidate_id=?',(cid,)).fetchone()
        if not force and old==(digest,VERSION):
            continue
        try: fields=json.loads(payload or '{}')
        except (TypeError,ValueError): fields={}
        r=classify(sentence,pep,sec,title)
        db.execute('''INSERT INTO processed_evidence VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
            ON CONFLICT(candidate_id) DO UPDATE SET
            paper_key=excluded.paper_key,peptide=excluded.peptide,function_code=excluded.function_code,
            source_section=excluded.source_section,doi=excluded.doi,paper_url=excluded.paper_url,
            source_text_hash=excluded.source_text_hash,processor_version=excluded.processor_version,
            result_json=excluded.result_json,extraction_status=excluded.extraction_status,processed_at=excluded.processed_at''',
            (cid,key,pep,clean(fields.get('Target_Function_Code')),sec,norm_doi(doi),url,
             digest,VERSION,json.dumps(r,ensure_ascii=False),r['extraction_status'],datetime.now(timezone.utc).isoformat()))
        count+=1
    db.commit()
    return count


def report(db):
    n=db.execute('SELECT COUNT(*) FROM processed_evidence').fetchone()[0]
    print('processed_evidence',n)
    for status,cnt in db.execute('SELECT extraction_status,COUNT(*) FROM processed_evidence GROUP BY extraction_status ORDER BY extraction_status'):
        print(status,cnt)
    reasons=Counter()
    for (raw,) in db.execute('SELECT result_json FROM processed_evidence WHERE extraction_status="REVIEW_REQUIRED"'):
        try: reasons[json.loads(raw).get('review_reason','미기재')]+=1
        except (ValueError,TypeError): pass
    print('=== 검토 사유 TOP 8 ===')
    for why,cnt in reasons.most_common(8): print(cnt,why)
    print('주의: 모든 자동 추출 결과는 원저 검증 전이며 추천 점수에 자동 반영되지 않습니다.')


def self_test():
    examples=[
      ('GHK increased wound closure in the treatment group to 60% compared with the control group at 40%.','GHK','','AUTO_EXTRACTED_PAIRED_VALUES',50.0),
      ('GHK increased wound closure by 30% compared with control.','GHK','','AUTO_EXTRACTED_REPORTED_CHANGE',None),
      ('GHK enhanced cell migration.','GHK','','QUALITATIVE_QUOTE',None),
      ('GHK did not improve wound closure.','GHK','','REVIEW_REQUIRED',None),
      ('GHK not only enhances biosafety but also promotes collagen synthesis.','GHK','','QUALITATIVE_QUOTE',None),
      ('RADA16-RGD promoted axonal regeneration [172].','RADA16-I','','REVIEW_REQUIRED',None),
      ('A different hydrogel increased cell migration.','GHK','','REVIEW_REQUIRED',None),
      ('RADA16-IKVAV promoted Schwann cell infiltration.','RADA16-I','','QUALITATIVE_QUOTE',None),
      ('GHK enhanced wound closure.','GHK','A review of regenerative peptides','REVIEW_REQUIRED',None),
    ]
    for s,p,t,expected,amount in examples:
        got=classify(s,p,title=t)
        assert got['extraction_status']==expected,(s,got)
        assert got['relative_change_percent']==amount,(s,got)
        assert not got['verified_for_ranking']
    print(f'SELF-TEST PASS ({len(examples)} synthetic sentences; original articles not verified)')


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--db',type=Path,default=DEFAULT_DB)
    parser.add_argument('--self-test',action='store_true')
    parser.add_argument('--report',action='store_true')
    parser.add_argument('--force',action='store_true')
    args=parser.parse_args()
    if args.self_test:
        self_test(); return
    if not args.db.is_file():
        parser.error('DB 파일을 찾을 수 없습니다: '+str(args.db))
    with sqlite3.connect(args.db) as db:
        initialize(db)
        if not args.report:
            print('NEW_OR_UPDATED',process(db,args.force))
        report(db)

if __name__=='__main__': main()

