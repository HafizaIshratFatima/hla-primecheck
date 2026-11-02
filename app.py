import re
from pathlib import Path
import pandas as pd
import streamlit as st

st.set_page_config(page_title="HLA PrimeCheck", page_icon="🧬", layout="wide")
st.title("🧬 HLA PrimeCheck")
st.caption("Allele-resolved sequence compatibility · IGA-050 and IGA-051 · Research prototype, not an editing-efficiency predictor")
ROOT = Path(__file__).parent

@st.cache_data
def load_data():
    ref = pd.read_csv(ROOT / "HLA_pegRNA_allele_specific_reference.csv", dtype={"Allele": str}).fillna("")
    exp = pd.read_csv(ROOT / "HLA_pegRNA_experimental_validation.csv", dtype=str).fillna("")
    return ref, exp

try:
    df, experiments = load_data()
except Exception as exc:
    st.error(f"Cannot load the reference data: {exc}")
    st.stop()

st.write("Enter **one or two alleles per locus**, separated by commas. Two-field typing (e.g. `A*24:02`) is supported, but may map to multiple higher-resolution database alleles.")
col_a, col_b, col_c = st.columns(3)
with col_a: inp_a = st.text_input("HLA-A", "A*24:02")
with col_b: inp_b = st.text_input("HLA-B", "B*52:01")
with col_c: inp_c = st.text_input("HLA-C", "C*12:02")

def normalize(value, locus):
    value = value.strip().upper().replace("HLA-", "").replace(" ", "")
    if not value:
        return ""
    if value.startswith("*"):
        value = locus + value
    elif value[0].isdigit():
        value = locus + "*" + value
    return value

def lookup(allele):
    exact = df[df["Allele"].str.upper() == allele].copy()
    if not exact.empty:
        return exact, "Exact database allele name"
    # For 2-field queries, keep higher-resolution records separate until summarizing.
    if re.fullmatch(r"[ABC]\*\d{2,3}:\d{2,3}", allele):
        sub = df[df["Allele"].str.upper().str.startswith(allele + ":")].copy()
        if not sub.empty:
            return sub, "Higher-resolution suballeles (inferred lookup)"
    return df.iloc[:0].copy(), "Not found in screened reference"

def category(row, i):
    prefix = f"pegRNA_{i}_"
    try:
        spacer = int(row[prefix + "spacer_mm"])
        pbs = int(row[prefix + "PBS_mm"])
        rtt = int(row[prefix + "RTT_mm"])
        pam = str(row[prefix + "PAM"])
    except (TypeError, ValueError):
        return "Incomplete sequence data"
    if spacer == 0 and pam == "NGG" and pbs == 0 and rtt == 0:
        return "Exact sequence-compatible"
    if spacer == 0 and pam == "NGG" and pbs == 0:
        return "Exact spacer/PAM/PBS; RTT near-match"
    return "Candidate with spacer/PAM/PBS differences"

def mm_range(values):
    numeric = pd.to_numeric(pd.Series(values), errors="coerce").dropna()
    if numeric.empty:
        return "Not assessed"
    low, high = int(numeric.min()), int(numeric.max())
    return str(low) if low == high else f"{low}–{high}"

inputs = []
for locus, raw in [("A", inp_a), ("B", inp_b), ("C", inp_c)]:
    for entry in re.split(r"[,;\n]+", raw):
        if entry.strip():
            inputs.append((locus, normalize(entry, locus)))

if st.button("Check allele compatibility", type="primary"):
    summary_rows = []
    detail_rows = []
    experimental_rows = []
    for locus, allele in inputs:
        if not re.fullmatch(rf"{locus}\*\d{{2,3}}:\d{{2,3}}(?::\d{{2,3}})*(?:[A-Z])?", allele):
            st.warning(f"Please check the format of `{allele}` (expected HLA-{locus} allele).")
            continue
        matched, lookup_status = lookup(allele)
        st.subheader(f"HLA-{locus}: {allele}")
        if matched.empty:
            st.warning("No matching record in the screened reference. This does **not** establish resistance to editing or absence of off-target risk.")
            for i, peg in [(1, "IGA-050"), (2, "IGA-051")]:
                summary_rows.append({"Entered allele": allele, "pegRNA": peg, "Database entries": 0, "Compatibility summary": "Not assessed", "Spacer MM": "—", "PBS MM": "—", "RTT MM": "—", "PAM(s)": "—", "Lookup": lookup_status})
        else:
            st.caption(f"{len(matched):,} database allele record(s) · {lookup_status}")
            for i, peg in [(1, "IGA-050"), (2, "IGA-051")]:
                categories = [category(r, i) for _, r in matched.iterrows()]
                distinct = sorted(set(categories))
                if len(distinct) == 1:
                    label = distinct[0]
                else:
                    label = "Variable among suballeles — expand details"
                summary_rows.append({
                    "Entered allele": allele,
                    "pegRNA": peg,
                    "Database entries": len(matched),
                    "Compatibility summary": label,
                    "Spacer MM": mm_range(matched[f"pegRNA_{i}_spacer_mm"]),
                    "PBS MM": mm_range(matched[f"pegRNA_{i}_PBS_mm"]),
                    "RTT MM": mm_range(matched[f"pegRNA_{i}_RTT_mm"]),
                    "PAM(s)": ", ".join(sorted(set(matched[f"pegRNA_{i}_PAM"].astype(str)))),
                    "Lookup": lookup_status,
                })
                for _, r in matched.iterrows():
                    detail_rows.append({
                        "Entered allele": allele, "Database allele": r["Allele"],
                        "pegRNA": peg, "Spacer MM": r[f"pegRNA_{i}_spacer_mm"],
                        "PAM": r[f"pegRNA_{i}_PAM"], "PBS MM": r[f"pegRNA_{i}_PBS_mm"],
                        "RTT MM": r[f"pegRNA_{i}_RTT_mm"],
                        "RTT differences": r[f"pegRNA_{i}_RTT_differences"],
                        "Sequence classification": category(r, i), "Lookup": lookup_status,
                    })
            st.dataframe(pd.DataFrame(summary_rows[-2:]), use_container_width=True, hide_index=True)
            if len(distinct) > 1 or lookup_status.startswith("Higher"):
                st.info("Two-field typing can represent multiple higher-resolution alleles. The summary describes database records, not a confirmed sequence for your particular cell line.")
            with st.expander(f"View {allele} higher-resolution sequence details ({len(matched)} database alleles)"):
                st.dataframe(pd.DataFrame(detail_rows[-2*len(matched):]), use_container_width=True, hide_index=True)
        # Experimental results only for the explicitly tested two-field typing.
        reported = "HLA-" + allele
        ev = experiments[experiments["Allele (reported experimental typing)"].str.upper() == reported.upper()]
        if not ev.empty:
            for _, r in ev.iterrows():
                experimental_rows.append({"Entered allele": allele, "Cell line": r["Cell line"], "IGA-050 editing (%)": r["pegRNA 1 editing %"], "IGA-051 editing (%)": r["pegRNA 2 editing %"], "Evidence": "Measured in the specified cell line; not transferable to all suballeles"})

    if summary_rows:
        st.divider()
        st.subheader("Downloadable cell-line summary")
        summary_df = pd.DataFrame(summary_rows)
        st.dataframe(summary_df, use_container_width=True, hide_index=True)
        st.download_button("Download summary CSV", summary_df.to_csv(index=False).encode("utf-8"), "HLA_PrimeCheck_summary.csv", "text/csv")
    if experimental_rows:
        st.subheader("Experimental validation (separate from sequence predictions)")
        st.dataframe(pd.DataFrame(experimental_rows).drop_duplicates(), use_container_width=True, hide_index=True)
        st.caption("Observed percentages are specific to the cited experiments and assay context; they are not predictions for other cell lines.")
    if detail_rows:
        with st.expander("All allele-resolved results (download and inspect)"):
            detailed_df = pd.DataFrame(detail_rows)
            st.dataframe(detailed_df, use_container_width=True, hide_index=True)
            st.download_button("Download detailed CSV", detailed_df.to_csv(index=False).encode("utf-8"), "HLA_PrimeCheck_detailed.csv", "text/csv")
    st.warning("Sequence compatibility is not a measured editing probability. Near-matching HLA-C sites may still be edited; a database non-hit does not establish specificity or safety.")
else:
    st.info("Click **Check allele compatibility** to view a concise summary and expandable allele details.")

with st.expander("Methods, evidence, and limitations"):
    st.markdown("""- Source: allele-specific reference derived from IPD-IMGT/HLA FASTA screening; verify database release and reference provenance before publication.
- Candidates were collected with up to three spacer mismatches and NGG/NAG/NGA PAMs. This is **not** a genome-wide off-target analysis.
- PBS/RTT mismatches are sequence comparisons under an assumed SpCas9 H840A nick and a **+CTAG insertion immediately after the nick**; validate against the pegRNA designs.
- For two-field input, matching higher-resolution alleles are grouped for display, **not assumed identical**. Mixed classifications are flagged.
- Experimental editing percentages come from specific tested alleles and cell lines; do not generalize them to every suballele or experimental condition.
- Absence of an allele in this screened table does not imply no editing or no off-target activity.
""")
st.caption("For research use only · Independently validate sequences, editing activity, and genomic integrity before publication or use.")
