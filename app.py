import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import math
import tempfile
import os
from paternity_logic import process_dataframe, parse_vcf
from str_logic import process_str_dataframe, parse_genemapper_export

st.set_page_config(
    page_title="Paternity & DNA Analysis Dashboard",
    page_icon="🧬",
    layout="wide"
)

# Custom Styling
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap');
    
    .main {
        background-color: #f8f9fc;
    }
    .stApp {
        font-family: 'Inter', sans-serif;
    }
    .metric-card {
        background: white;
        border-radius: 12px;
        padding: 24px;
        box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05), 0 2px 4px -1px rgba(0,0,0,0.03);
        text-align: center;
        border: 1px solid #e2e8f0;
        transition: transform 0.2s ease, box-shadow 0.2s ease;
    }
    .metric-card:hover {
        transform: translateY(-2px);
        box-shadow: 0 10px 15px -3px rgba(0,0,0,0.1), 0 4px 6px -2px rgba(0,0,0,0.05);
    }
    .metric-title {
        color: #64748b;
        font-size: 13px;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.05em;
    }
    .metric-value {
        color: #0f172a;
        font-size: 36px;
        font-weight: 700;
        margin-top: 8px;
        background: -webkit-linear-gradient(45deg, #2563eb, #3b82f6);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
    .metric-value-neutral {
        color: #0f172a;
        font-size: 36px;
        font-weight: 700;
        margin-top: 8px;
    }
    .result-text {
        font-size: 24px;
        font-weight: 700;
        text-align: center;
        padding: 20px;
        border-radius: 12px;
        margin-top: 24px;
        box-shadow: 0 4px 6px -1px rgba(0,0,0,0.05);
    }
    .result-strong {
        background: linear-gradient(to right, #dcfce7, #bbf7d0);
        color: #166534;
        border: 1px solid #86efac;
    }
    .result-weak {
        background: linear-gradient(to right, #fef9c3, #fef08a);
        color: #854d0e;
        border: 1px solid #fde047;
    }
    .result-non {
        background: linear-gradient(to right, #fee2e2, #fecaca);
        color: #991b1b;
        border: 1px solid #fca5a5;
    }
    .app-title {
        font-weight: 800;
        background: -webkit-linear-gradient(45deg, #1e3a8a, #3b82f6);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
    }
</style>
""", unsafe_allow_html=True)

st.markdown('<h1 class="app-title">🧬 Paternity & DNA Analysis Dashboard</h1>', unsafe_allow_html=True)
st.markdown("Platform komprehensif untuk uji paternitas, mendukung analisis genetik melalui marker *Short Tandem Repeat* (STR) maupun *Single Nucleotide Polymorphism* (SNP).")
st.markdown("<br>", unsafe_allow_html=True)

def get_interpretation_class(interp_text):
    if "NON-PATERNITAS" in interp_text or "NON-paternitas" in interp_text or "Tidak ada" in interp_text:
        return "result-non"
    elif "SANGAT KUAT" in interp_text or "KUAT" in interp_text or "TERBUKTI" in interp_text:
        return "result-strong"
    else:
        return "result-weak"

tab1, tab2 = st.tabs(["🔬 Uji Short Tandem Repeat (Standar Medis)", "🧬 Uji Profil Lanjutan (SNP)"])

with tab1:
    st.header("Analisis Paternitas Konklusif (STR)")
    st.markdown("STR adalah standar emas di dunia hukum dan medis untuk pengujian genetik. Unggah file CSV (*format:* `Locus`, `Mother`, `Alleged_Father`, `Child`) atau **Ekspor GeneMapper (.txt)**.")
    
    uploaded_file_str = st.file_uploader("📂 Unggah Dataset STR", type=["csv", "txt"], key="file_str")
    
    if uploaded_file_str is not None:
        try:
            # Handle file naming and parsing
            file_details = {"FileName": uploaded_file_str.name, "FileType": uploaded_file_str.type}
            
            if uploaded_file_str.name.endswith('.txt'):
                # Save to a secure temporary file
                with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as tmp:
                    tmp.write(uploaded_file_str.getbuffer())
                    tmp_path = tmp.name
                try:
                    df_str = parse_genemapper_export(tmp_path)
                finally:
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
            else:
                df_str = pd.read_csv(uploaded_file_str)
            required_cols = ["Locus", "Mother", "Alleged_Father", "Child"]
            
            if not all(col in df_str.columns for col in required_cols):
                st.error(f"Format CSV STR tidak sesuai! Kolom yang dibutuhkan: {', '.join(required_cols)}")
            else:
                with st.spinner('Menghitung Paternity Index lokus STR...'):
                    cpi, prob_paternity, total_loci, matching_loci, ex_loci, interp_str, layman_str, df_res = process_str_dataframe(df_str)
                
                if total_loci > 0:
                    col1, col2, col3, col4 = st.columns(4)
                    with col1:
                        st.markdown(f'<div class="metric-card"><div class="metric-title">Total Lokus</div><div class="metric-value-neutral">{total_loci}</div></div>', unsafe_allow_html=True)
                    with col2:
                        st.markdown(f'<div class="metric-card"><div class="metric-title">Lokus Identik</div><div class="metric-value" style="background: -webkit-linear-gradient(45deg, #10b981, #34d399); -webkit-background-clip: text;">{matching_loci}</div></div>', unsafe_allow_html=True)
                    with col3:
                        st.markdown(f'<div class="metric-card"><div class="metric-title">Lokus Eksklusi</div><div class="metric-value" style="background: -webkit-linear-gradient(45deg, #ef4444, #f87171); -webkit-background-clip: text;">{ex_loci}</div></div>', unsafe_allow_html=True)
                    with col4:
                        st.markdown(f'<div class="metric-card"><div class="metric-title">Probabilitas (W)</div><div class="metric-value">{prob_paternity:.4f}%</div></div>', unsafe_allow_html=True)
                        
                    st.markdown(f'<div class="result-text {get_interpretation_class(interp_str)}">🎯 KESIMPULAN: {interp_str}</div>', unsafe_allow_html=True)
                    st.info(f"**📚 Penjelasan Medis (Awam):** {layman_str}")
                    
                    st.markdown("---")
                    st.subheader("Data Rinci Per Lokus STR")
                    
                    def style_str_dataframe(row):
                        if "Eksklusi" in str(row['Status STR']) or "Exclusion" in str(row['Status STR']):
                            return ['background-color: #fee2e2; color: #991b1b'] * len(row)
                        elif "Inklusi" in str(row['Status STR']):
                            return ['background-color: #dcfce7; color: #166534'] * len(row)
                        return [''] * len(row)
                        
                    st.dataframe(df_res.style.apply(style_str_dataframe, axis=1), use_container_width=True, height=350)
                    
                    st.markdown(f"*Keterangan Khusus:* Cumulative Paternity Index (CPI) Anda bernilai eksak **{cpi:,.2f}**.")
                    
                    csv_export_str = df_res.to_csv(index=False).encode('utf-8')
                    st.download_button(label="📥 Unduh Laporan Medis STR (CSV)", data=csv_export_str, file_name='Laporan_Medis_STR.csv', mime='text/csv')
                        
        except Exception as e:
            st.error(f"Terjadi kesalahan teknis pemrosesan: {str(e)}")
    else:
        st.info("💡 **Tips:** Anda dapat mencoba menggunakan file `str_sample.csv` yang baru saja sistem generate untuk melihat simulasi hasil konklusif!")

with tab2:
    st.header("Analisis Ekstensif Profil SNP")
    st.markdown("Metode eksplorasi pendukung menggunakan Single Nucleotide Polymorphism. Pastikan format kolom: `SNP_ID`, `Mother`, `Alleged_Father`, `Child` atau gunakan **Format Industri VCF (.vcf)**.")
    
    uploaded_file_snp = st.file_uploader("📂 Unggah Dataset SNP", type=["csv", "vcf"], key="file_snp")
    
    if uploaded_file_snp is not None:
        try:
            if uploaded_file_snp.name.endswith('.vcf'):
                # Save to a secure temporary file
                with tempfile.NamedTemporaryFile(delete=False, suffix=".vcf") as tmp:
                    tmp.write(uploaded_file_snp.getbuffer())
                    tmp_path = tmp.name
                try:
                    df_snp = parse_vcf(tmp_path)
                finally:
                    if os.path.exists(tmp_path):
                        os.remove(tmp_path)
            else:
                df_snp = pd.read_csv(uploaded_file_snp)
            
            required_cols = ["SNP_ID", "Mother", "Alleged_Father", "Child"]
            if not all(col in df_snp.columns for col in required_cols):
                st.error(f"Format CSV SNP tidak sesuai! Kolom yang dibutuhkan: {', '.join(required_cols)}")
            else:
                with st.spinner('Menghitung kesimpulan likelihood (HWE)...'):
                    gen_ibu, gen_ayah, gen_anak, snps_tidak_hwe, L_pat, L_non, LR, interpretasi, penjelasan_awam = process_dataframe(df_snp)
                
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.markdown(f'<div class="metric-card"><div class="metric-title">Total Cek SNP</div><div class="metric-value-neutral">{len(df_snp)}</div></div>', unsafe_allow_html=True)
                with col2:
                    valid_count = len(df_snp) - len(snps_tidak_hwe)
                    st.markdown(f'<div class="metric-card"><div class="metric-title">SNP Valid (HWE)</div><div class="metric-value-neutral">{valid_count}</div></div>', unsafe_allow_html=True)
                with col3:
                    st.markdown(f'<div class="metric-card"><div class="metric-title">Likelihood Ratio</div><div class="metric-value">{LR:.2e}</div></div>', unsafe_allow_html=True)
                    
                st.markdown(f'<div class="result-text {get_interpretation_class(interpretasi)}">🔍 {interpretasi}</div>', unsafe_allow_html=True)
                st.info(f"**📚 Analisis Risiko:** {penjelasan_awam}")
                
                with st.expander("📊 Tampilkan Rincian Data Panel & Grafik HWE"):
                    col_chart1, col_chart2 = st.columns(2)
                    with col_chart1:
                        pie_data = pd.DataFrame({
                            "Status": ["Valid (Memenuhi HWE)", "Tidak Valid (Gagal HWE)"],
                            "Jumlah": [valid_count, len(snps_tidak_hwe)]
                        })
                        fig_pie = px.pie(pie_data, values="Jumlah", names="Status", color="Status", color_discrete_map={"Valid (Memenuhi HWE)": "#10b981", "Tidak Valid (Gagal HWE)": "#ef4444"}, hole=0.4)
                        st.plotly_chart(fig_pie, use_container_width=True)
                        
                    with col_chart2:
                        try: log_lr = math.log10(LR) if LR > 0 and LR != float('inf') else 0
                        except: log_lr = 0
                        fig_gauge = go.Figure(go.Indicator(
                            mode = "gauge+number", value = log_lr, title = {'text': "Log10(Ratio)"},
                            gauge = {'axis': {'range': [-2, 5]}, 'bar': {'color': "#3b82f6"}, 'steps': [
                                    {'range': [-2, 0], 'color': "#fee2e2"}, {'range': [2, 5], 'color': "#22c55e"}]}))
                        st.plotly_chart(fig_gauge, use_container_width=True)
                        
                    df_snp['Status HWE'] = df_snp['SNP_ID'].apply(lambda x: "❌ Gagal" if str(x) in snps_tidak_hwe else "✅ Valid")
                    st.dataframe(df_snp, height=250)
                    
        except Exception as e:
            st.error(f"Terjadi kesalahan teknis: {str(e)}")
    else:
        st.info("💡 Unggah dataset SNP Anda di sini jika uji STR dirasa kurang bervariasi atau data sampel rusak secara biomolekuler (degradasi DNA).")
