"""Tissue Regeneration Peptide Navigator — conservative Europe PMC builder.

Original Excel is READ ONLY. Automated sentences are NEVER verified evidence.
Run: python build_database.py --peptide GHK --limit 10
     python build_database.py --all --limit 10
     python build_database.py --report
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import sqlite3
import time
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

import pandas as pd

ROOT = Path(__file__).resolve().parent
EXCEL = ROOT / 'peptide_regeneration_DB_3차_screening.xlsx'
DB = ROOT / 'navigator.sqlite3'
API = 'https://www.ebi.ac.uk/europepmc/webservices/rest'
FUNCTIONS = {
    'WOUND': ('wound', 'healing', 'epithelial', 'scratch', 'wound closure'),
    'ANGIO': ('angiogenesis', 'endothelial', 'tube formation', 'vascular'),
    'NERVE': ('nerve regeneration', 'neurite', 'axon', 'schwann', 'neural'),
}
ALIASES = {
    'GHK': ['"GHK"', '"glycyl-histidyl-lysine"', '"GHK-Cu"'],
    'RADA16-I': ['"RADA16"', '"RADA16-I"'],
    'QK': ['"QK peptide"', '"VEGF mimetic peptide"'],
    'IKVAV': ['"IKVAV"'],
    'KLD-12': ['"KLD12"', '"KLD-12"'],
    'RGD': ['"RGD peptide"', '"arginine-glycine-aspartate"'],
    'YIGSR': ['"YIGSR"'],
}
NA = 'Not reported'


def norm(value):
    return re.sub(r'\s+', ' ', str(value or '')).strip()


def doi_norm(value):
    return re.sub(r'^https?://(?:dx\.)?doi\.org/', '', norm(value), flags=re.I).lower()


def now():
    return datetime.now(timezone.utc).isoformat()


def get(url, params=None, as_xml=False):
    if params:
        url += '?' + urlencode(params)
    req = Request(url, headers={'User-Agent': 'TissueRegenerationNavigator/0.1 (academic research)', 'Accept': 'application/xml' if as_xml else 'application/json'})
    for attempt in range(3):
        try:
            with urlopen(req, timeout=25) as resp:
                data = resp.read(10_000_001)
            if len(data) > 10_000_000:
                raise ValueError('API result too large')
            return ET.fromstring(data) if as_xml else json.loads(data)
        except HTTPError as exc:
            if exc.code in (400, 401, 403, 404) or attempt == 2:
                raise
        except (URLError, TimeoutError):
            if attempt == 2:
                raise
        time.sleep(2 ** attempt)


def initialize(db):
    db.executescript('''
    PRAGMA foreign_keys=ON;
    CREATE TABLE IF NOT EXISTS excel_screening (
      record_id TEXT PRIMARY KEY, peptide TEXT, function_code TEXT,
      study_id TEXT, doi TEXT, payload_json TEXT NOT NULL
    );
    CREATE TABLE IF NOT EXISTS excel_recommendations (
      recommendation_id TEXT PRIMARY KEY, peptide_id TEXT, function_code TEXT,
      doi TEXT, payload_json TEXT NOT NULL
    );
    CREATE TABLE IF NOT EXISTS papers (
      paper_key TEXT PRIMARY KEY, title TEXT, year TEXT, doi TEXT,
      pmid TEXT, pmcid TEXT, paper_url TEXT, abstract TEXT,
      fulltext_status TEXT, collected_at TEXT
    );
    CREATE TABLE IF NOT EXISTS paper_peptides (
      paper_key TEXT NOT NULL REFERENCES papers(paper_key),
      peptide TEXT NOT NULL,
      PRIMARY KEY(paper_key,peptide)
    );
    CREATE TABLE IF NOT EXISTS extraction_candidates (
      candidate_id TEXT PRIMARY KEY, paper_key TEXT NOT NULL REFERENCES papers(paper_key),
      peptide TEXT, source_section TEXT, source_text TEXT NOT NULL,
      screening_json TEXT NOT NULL, evidence_status TEXT NOT NULL,
      created_at TEXT NOT NULL
    );
    CREATE INDEX IF NOT EXISTS idx_candidate_peptide ON extraction_candidates(peptide);
    ''')
    db.commit()


def import_excel(db, path):
    required = ['Screening DB v3', 'Recommendation DB v3', 'Peptide Library']
    sheets = pd.read_excel(path, sheet_name=required, dtype=str, engine='openpyxl', keep_default_na=False)
    library = sheets['Peptide Library']
    available = {norm(x) for x in library['Name']}
    if not set(ALIASES).issubset(available):
        raise ValueError(f'Peptide Library mismatch: {sorted(set(ALIASES) - available)}')
    screening = sheets['Screening DB v3']
    recommendation = sheets['Recommendation DB v3']
    expected = {'Record_ID', 'Original_Peptide', 'Target_Function_Code', 'Study_ID', 'DOI'}
    if not expected.issubset(screening.columns):
        raise ValueError('Screening DB v3 has missing required columns')
    with db:
        for _, row in screening.iterrows():
            item = row.to_dict()
            db.execute('''INSERT OR IGNORE INTO excel_screening
                (record_id,peptide,function_code,study_id,doi,payload_json)
                VALUES (?,?,?,?,?,?)''', (norm(item['Record_ID']),norm(item['Original_Peptide']),norm(item['Target_Function_Code']),norm(item['Study_ID']),doi_norm(item['DOI']),json.dumps(item,ensure_ascii=False)))
        for _, row in recommendation.iterrows():
            item = row.to_dict()
            db.execute('''INSERT OR IGNORE INTO excel_recommendations
                (recommendation_id,peptide_id,function_code,doi,payload_json)
                VALUES (?,?,?,?,?)''', (norm(item['Recommendation_ID']),norm(item['Peptide_ID']),norm(item['Target_Function_Code']),doi_norm(item['DOI']),json.dumps(item,ensure_ascii=False)))
    return len(screening), len(recommendation)


def paper_key(p):
    doi = doi_norm(p.get('doi'))
    if doi: return 'DOI:' + doi
    if p.get('pmid'): return 'PMID:' + norm(p['pmid'])
    if p.get('pmcid'): return 'PMCID:' + norm(p['pmcid'])
    return ''  # Do not invent a title-based scientific identifier.


def query_for(peptide):
    return '(' + ' OR '.join(ALIASES[peptide]) + ') AND (wound OR regenerat* OR angiogen* OR nerve OR tissue OR hydrogel) AND (SRC:MED OR SRC:PMC)'


def sentences(text):
    # The segmentation may be imperfect; all results remain UNVERIFIED.
    return [norm(s) for s in re.split(r'(?<=[.!?])\s+(?=[A-Z(])', norm(text)) if 35 <= len(norm(s)) <= 1200]


def candidates_from_text(text, peptide, section):
    out = []
    # DO NOT infer intervention-specific efficacy merely from terms in a sentence.
    for sentence in sentences(text):
        if not re.search(r'\b(increas\w*|decreas\w*|enhanc\w*|reduc\w*|promot\w*|improv\w*|significan\w*|compar\w*)\b', sentence, re.I):
            continue
        hits = [code for code, terms in FUNCTIONS.items() if any(term in sentence.lower() for term in terms)]
        if not hits:
            continue
        for code in hits:
            item = {
                'Record_ID': NA, 'Original_Peptide': peptide,
                'Original_Sequence': NA, 'Modification_Position': NA,
                'Modification_Method': NA, 'Modified_Peptide': NA,
                'Modified_Sequence': NA, 'Formulation': NA,
                'Intervention_Class': '확인 필요', 'Target_Function': NA,
                'Target_Function_Code': code, 'Evaluation_Method': NA,
                'Measured_Outcome': NA, 'Effect_Direction': '확인 필요',
                'Quantitative_Result': NA, 'Statistical_Significance': NA,
                'Comparator': NA, 'Comparator_Type': NA,
                'Evidence_Status': 'UNVERIFIED_EXTRACTION',
                'Source_Verification': section,
                'Reference': NA, 'DOI': NA, 'Paper_URL': NA,
                'Study_ID': NA,
                'Candidate_Statement': sentence,
                'Reviewer_Notes': 'Keyword-selected sentence only. No assignment of intervention, comparator, direction or quantitative effect.'
            }
            out.append((sentence,code,item))
    return out


def save_paper(db, p, peptide, use_fulltext=True):
    key = paper_key(p)
    if not key: return 0, 0, 'missing_paper_id'
    pmcid = norm(p.get('pmcid'))
    doi = doi_norm(p.get('doi'))
    abstract = norm(p.get('abstractText'))
    url = 'https://europepmc.org/articles/' + pmcid if pmcid else 'https://doi.org/' + doi if doi else 'https://europepmc.org/article/MED/' + norm(p.get('pmid'))
    exists = db.execute('SELECT paper_key FROM papers WHERE paper_key=?',(key,)).fetchone()
    if exists:
        with db:
            db.execute('INSERT OR IGNORE INTO paper_peptides VALUES (?,?)',(key,peptide))
        return 0,0,'cached'
    sections = [('abstract',abstract)] if abstract else []
    fulltext_status = 'not_attempted'
    if use_fulltext and re.fullmatch(r'PMC\d+',pmcid):
        try:
            root = get(f'{API}/{pmcid}/fullTextXML',as_xml=True)
            sections.extend(('fulltext',norm(' '.join(el.itertext()))) for el in root.findall('.//body//sec/p'))
            fulltext_status = 'retrieved'
        except HTTPError as exc:
            fulltext_status = f'unavailable_http_{exc.code}'
        except (URLError, TimeoutError, ValueError, ET.ParseError) as exc:
            fulltext_status = 'fetch_error_' + type(exc).__name__
    else:
        fulltext_status = 'no_pmcid_or_disabled'
    count = 0
    with db:
        db.execute('''INSERT OR IGNORE INTO papers VALUES (?,?,?,?,?,?,?,?,?,?)''',
            (key,norm(p.get('title')),norm(p.get('pubYear')),doi,norm(p.get('pmid')),pmcid,url,abstract,fulltext_status,now()))
        db.execute('INSERT OR IGNORE INTO paper_peptides VALUES (?,?)',(key,peptide))
        for section,text in sections:
            for sentence,code,payload in candidates_from_text(text,peptide,section):
                payload.update({'Reference':norm(p.get('title')),'DOI':doi or NA,'Paper_URL':url,'Study_Year':norm(p.get('pubYear'))})
                seed = '\x1f'.join([key,peptide,code,section,sentence])
                cid = hashlib.sha256(seed.encode()).hexdigest()[:24]
                db.execute('''INSERT OR IGNORE INTO extraction_candidates VALUES (?,?,?,?,?,?,?,?)''',
                    (cid,key,peptide,section,sentence,json.dumps(payload,ensure_ascii=False),'UNVERIFIED_EXTRACTION',now()))
                count += 1
    return 1,count,fulltext_status


def collect(db, peptide, limit, fulltext=True):
    query = query_for(peptide)
    print('QUERY:', query)
    data = get(API + '/search', {'query':query,'format':'json','resultType':'core','pageSize':min(limit,100),'cursorMark':'*'})
    results = data.get('resultList',{}).get('result',[])
    created = candidates = 0
    stats = {}
    for paper in results[:limit]:
        n,c,status=save_paper(db,paper,peptide,fulltext)
        created += n; candidates += c
        stats[status] = stats.get(status,0)+1
        time.sleep(0.15)
    print(json.dumps({'peptide':peptide,'api_matches':data.get('hitCount',NA),'fetched':len(results[:limit]),'new_papers':created,'new_candidate_sentences_attempted':candidates,'statuses':stats},ensure_ascii=False))


def report(db):
    for table in ('excel_screening','excel_recommendations','papers','extraction_candidates'):
        print(table,db.execute('SELECT count(*) FROM '+table).fetchone()[0])
    print('UNVERIFIED:',db.execute("SELECT count(*) FROM extraction_candidates WHERE evidence_status='UNVERIFIED_EXTRACTION'").fetchone()[0])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--peptide',choices=sorted(ALIASES),default='GHK')
    parser.add_argument('--all',action='store_true')
    parser.add_argument('--limit',type=int,default=10,help='Max API papers per peptide; maximum 50')
    parser.add_argument('--no-fulltext',action='store_true')
    parser.add_argument('--report',action='store_true')
    parser.add_argument('--excel',type=Path,default=EXCEL)
    parser.add_argument('--db',type=Path,default=DB)
    args=parser.parse_args()
    if not 1 <= args.limit <= 50: parser.error('--limit must be 1..50')
    if not args.excel.is_file(): parser.error(f'Excel not found: {args.excel}')
    db=sqlite3.connect(args.db,timeout=30)
    try:
        initialize(db)
        print('IMPORTED',import_excel(db,args.excel))
        if not args.report:
            for peptide in (ALIASES if args.all else [args.peptide]):
                collect(db,peptide,args.limit,not args.no_fulltext)
        report(db)
    finally:db.close()

if __name__ == '__main__': main()
