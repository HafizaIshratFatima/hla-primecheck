# HLA PrimeCheck — research prototype

## Local test
Place `app.py`, `requirements.txt` and `HLA_pegRNA_allele_specific_reference.csv` in one folder. Then run:

```bash
python3 -m pip install -r requirements.txt
python3 -m streamlit run app.py
```

## Publish publicly using Streamlit Community Cloud
1. Create a GitHub account and a new repository, e.g. `hla-primecheck`.
2. Upload `app.py`, `requirements.txt`, and `HLA_pegRNA_allele_specific_reference.csv` to the repository root.
3. Go to https://share.streamlit.io/ and sign in with GitHub.
4. Choose **Create app**, select your repository, branch `main`, and main file `app.py`.
5. Deploy, then share the public Streamlit app URL.

Review the database, allele nomenclature, license/attribution, source release and methodology before publishing. A public GitHub repository exposes the data and code. This is an exploratory sequence lookup, not a clinical tool or editing-efficiency predictor.
