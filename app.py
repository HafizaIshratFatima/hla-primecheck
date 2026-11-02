import re
from pathlib import Path
import pandas as pd
import streamlit as st

st.set_page_config(page_title='HLA PrimeCheck', page_icon='🧬', layout='wide')
st.title('🧬 HLA PrimeCheck')
st.caption('Research prototype · IGA-050 (pegRNA 1) and IGA-051 (pegRNA 2) · sequence compatibility, not predicted editing efficiency')
DATA = Path(__file__).parent / 'HLA_pegRNA_allele_specific_reference.csv'
@st.cache_data
def load():
    return pd.read_csv(DATA, dtype={'Allele':str}).fillna('')
df = load()
st.markdown('Enter **one or two alleles per locus**, separated by commas. Use the full allele names where possible (for example, `A*24:02`).')
a,b,c = st.columns(3)
with a: aa=st.text_input('HLA-A', 'A*24:02')
with b: bb=st.text_input('HLA-B', 'B*52:01')
with c: cc=st.text_input('HLA-C', 'C*12:02')

def norm(v, locus):
    v=v.strip().upper().replace('HLA-','').replace(' ','')
    if not v: return ''
    if not v.startswith(locus+'*'):
        if v.startswith('*'): v=locus+v
        elif v[0].isdigit(): v=locus+'*'+v
    return v

def lookup(allele):
    exact=df[df['Allele'].str.upper()==allele]
    if len(exact): return exact.copy(), 'Exact allele-name match'
    # Two-field input: do not imply suballeles have identical target sequences.
    if re.fullmatch(r'[ABC]\*\d{2,3}:\d{2,3}',allele):
        sub=df[df['Allele'].str.upper().str.startswith(allele+':')].copy()
        if len(sub): return sub, 'Higher-resolution suballeles (NOT an exact two-field match)'
    return df.iloc[0:0].copy(), 'No matching entry in this screened database'

requests=[]
for locus,field in [('A',aa),('B',bb),('C',cc)]:
    requests.extend((locus,norm(x,locus)) for x in re.split(r'[,;\n]+',field) if x.strip())
if st.button('Check allele compatibility', type='primary'):
    results=[]
    for locus,allele in requests:
        matched,status=lookup(allele)
        if matched.empty:
            results.append({'Entered allele':allele,'Database allele':'—','Lookup':status,'pegRNA':'—','Spacer MM':'—','PAM':'—','PBS MM':'—','RTT MM':'—','Sequence classification':'Not assessed'})
            continue
        for _,r in matched.iterrows():
            for i,name in [(1,'IGA-050'),(2,'IGA-051')]:
                results.append({'Entered allele':allele,'Database allele':r['Allele'],'Lookup':status,'pegRNA':name,'Spacer MM':r[f'pegRNA_{i}_spacer_mm'],'PAM':r[f'pegRNA_{i}_PAM'],'PBS MM':r[f'pegRNA_{i}_PBS_mm'],'RTT MM':r[f'pegRNA_{i}_RTT_mm'],'Sequence classification':r[f'pegRNA_{i}_class']})
    out=pd.DataFrame(results)
    st.dataframe(out,use_container_width=True,hide_index=True)
    st.download_button('Download results as CSV',out.to_csv(index=False).encode(),'HLA_PrimeCheck_lookup.csv','text/csv')
    st.info('An exact sequence match is not proof of editing. Mismatched sites may still edit; absence from this screening does not establish safety. Two-field HLA types can encompass different higher-resolution sequences.')
else:
    st.write('Select **Check allele compatibility** to display results.')
with st.expander('Methods and limitations'):
    st.markdown('''- Source: allele-specific reference CSV derived from the supplied IPD-IMGT/HLA FASTA screening. Check release/version and sequence completeness before publication.
- Screened candidate sites have up to three spacer mismatches and recognized PAM classes (NGG, NAG, NGA). The tool is not a genome-wide off-target search.
- RTT mismatch values describe a **fixed +CTAG insertion immediately after the assumed nick**, and are sequence comparisons only.
- **Exact two-field allele inputs** are not automatically equated with higher-resolution suballeles. When no exact entry exists, suballeles are shown individually.
- Sequence classes and allele counts are not editing probabilities, population coverage, or clinical suitability.
- Experimental validation for A*24:02 and B*52:01 should be displayed separately with assay conditions, not propagated to all suballeles.''')
st.caption('For research use only · Verify candidate sequences and assay specificity experimentally.')
