import pandas as pd
from typing import Dict, Tuple, List, Any

# Konstanta frekuensi alel dasar populasi lokal (jika tidak disediakan di data)
DEFAULT_ALLELE_FREQ = 0.05

def parse_genemapper_export(file_path: str) -> pd.DataFrame:
    """
    Membaca file ekspor GeneMapper (.txt tab-delimited).
    Mengonversi format vertikal (per sampel) ke horizontal (trio).
    """
    try:
        # GeneMapper export biasanya tab-delimited
        df_raw = pd.read_csv(file_path, sep='\t')
        
        # Kolom standar GeneMapper: 'Sample Name', 'Marker', 'Allele 1', 'Allele 2'
        cols = df_raw.columns
        sample_col = 'Sample Name' if 'Sample Name' in cols else cols[0]
        marker_col = 'Marker' if 'Marker' in cols else cols[1]
        a1_col = 'Allele 1' if 'Allele 1' in cols else cols[2]
        a2_col = 'Allele 2' if 'Allele 2' in cols else cols[3]
        
        # Filter hanya untuk CHILD, MOTHER, FATHER
        target_samples = ['CHILD', 'MOTHER', 'FATHER']
        df_filtered = df_raw[df_raw[sample_col].isin(target_samples)].copy()
        
        # Gabungkan Allele 1 & 2 menjadi string 'A1, A2'
        df_filtered['Genotype'] = df_filtered[a1_col].astype(str) + ", " + df_filtered[a2_col].astype(str)
        
        # Pivot agar Marker menjadi baris, dan Sample Name menjadi kolom
        df_pivoted = df_filtered.pivot(index=marker_col, columns=sample_col, values='Genotype')
        
        # Reset index agar 'Marker' jadi kolom biasa (Locus)
        df_pivoted = df_pivoted.reset_index().rename(columns={marker_col: 'Locus'})
        
        # Pastikan kolom wajib ada
        for s in target_samples:
            if s not in df_pivoted.columns:
                df_pivoted[s] = ""
                
        # Rename FATHER -> Alleged_Father
        df_pivoted = df_pivoted.rename(columns={'FATHER': 'Alleged_Father'})
        
        return df_pivoted
    except Exception as e:
        print(f"Error parsing GeneMapper export: {e}")
        return pd.DataFrame()

def parse_alleles(allele_str: str) -> List[float]:
    """Parse string allele seperti '10, 11' atau '9.3/10' menjadi list float"""
    if pd.isna(allele_str):
        return []
    # Bersihkan string dan split berdasarkan koma atau slash
    cleaned = str(allele_str).replace(' ', '')
    if ',' in cleaned:
        parts = cleaned.split(',')
    elif '/' in cleaned:
        parts = cleaned.split('/')
    else:
        parts = [cleaned, cleaned] # Homozygous jika hanya 1 angka
    
    try:
        return [float(p) for p in parts]
    except ValueError:
        return []

def calculate_str_locus_pi(child_alleles: List[float], mother_alleles: List[float], father_alleles: List[float], pop_freq: float = DEFAULT_ALLELE_FREQ) -> Tuple[float, str]:
    """
    Hitung Paternity Index (PI) untuk satu lokus STR.
    Return: (PI, status_match)
    """
    if len(child_alleles) < 2:
        child_alleles = [child_alleles[0], child_alleles[0]] if child_alleles else []
    if len(mother_alleles) < 2:
        mother_alleles = [mother_alleles[0], mother_alleles[0]] if mother_alleles else []
    if len(father_alleles) < 2:
        father_alleles = [father_alleles[0], father_alleles[0]] if father_alleles else []
        
    if not child_alleles or not mother_alleles or not father_alleles:
         return 0.0, "⚠️ Data Tidak Lengkap"

    # Cari kandidat alel ayah (obligate paternal allele)
    c1, c2 = child_alleles[0], child_alleles[1]
    m1, m2 = mother_alleles[0], mother_alleles[1]
    
    paternal_candidates = []
    
    # Skenario 1: c1 dari ibu, c2 dari ayah
    if c1 in [m1, m2]:
        paternal_candidates.append(c2)
    # Skenario 2: c2 dari ibu, c1 dari ayah
    if c2 in [m1, m2]:
        paternal_candidates.append(c1)
        
    if not paternal_candidates:
         # Kemungkinan Mutasi pada Ibu atau Bukan Ibu Kandung
        return 0.001, "❌ Maternal Exclusion (Cek Mutasi Ibu)"
        
    # Periksa apakah terduga ayah memiliki alel kandidat
    match = False
    for candidate in paternal_candidates:
        if candidate in father_alleles:
            match = True
            break
            
    if not match:
        # Pria ini tidak punya alel yang diwajibkan anak.
        # Bisa jadi karena mutasi ayah, sehingga sering diberikan PI mutasi standar ~0.001
        return 0.001, "❌ Eksklusi"
        
    # Jika homozigot, sumbangannya = 1/freq. Jika heterozigot = 1/2freq
    if father_alleles[0] == father_alleles[1]:
        pi = 1.0 / pop_freq
    else:
        pi = 1.0 / (2 * pop_freq)
        
    return pi, "✅ Inklusi"

def process_str_dataframe(df: pd.DataFrame) -> Tuple[float, float, int, int, int, str, str, pd.DataFrame]:
    """
    Memproses DataFrame Lokus STR.
    Kolom Wajib: Locus, Mother, Alleged_Father, Child. (Allele_Freq bersifat opsional).
    """
    cpi = 1.0
    matching_loci = 0
    excluded_loci = 0
    
    results = []
    
    for _, row in df.iterrows():
        try:
            locus = str(row["Locus"]).strip()
        except KeyError:
            continue
            
        freq = row.get("Allele_Freq", DEFAULT_ALLELE_FREQ)
        
        m_alleles = parse_alleles(row.get("Mother", ""))
        f_alleles = parse_alleles(row.get("Alleged_Father", ""))
        c_alleles = parse_alleles(row.get("Child", ""))
        
        # Validasi format sebelum menghitung PI
        if not m_alleles or not f_alleles or not c_alleles:
            continue
            
        pi, status = calculate_str_locus_pi(c_alleles, m_alleles, f_alleles, float(freq) if pd.notna(freq) else DEFAULT_ALLELE_FREQ)
        
        if "Inklusi" in status:
            matching_loci += 1
            cpi *= pi
        elif "Eksklusi" in status:
            excluded_loci += 1
            cpi *= 0.001  # Mutasi index adjustment
            
        results.append({
            "Lokus (Marker)": locus,
            "Ibu": row.get("Mother", ""),
            "Terduga Ayah": row.get("Alleged_Father", ""),
            "Anak": row.get("Child", ""),
            "Status STR": status,
            "PI (Paternity Index)": round(pi, 4) if pi > 0.001 else "<0.001"
        })
        
    # Probability of Paternity (W) = CPI / (CPI + 1) -> Assuming prior probability of 50%
    prob_paternity = (cpi / (cpi + 1)) * 100 if cpi > 0 else 0.0
    total_loci = len(results)
    
    if total_loci == 0:
         return 0.0, 0.0, 0, 0, 0, "DATA TIDAK VALID", "Tidak dapat membaca baris marker. Pastikan format kolom sesuai dengan template.", pd.DataFrame()
         
    # Interpretasi STR di dunia nyata (berdasarkan Standar AABB)
    if excluded_loci >= 3:
        interp = "NON-PATERNITAS"
        layman = f"Sebanyak {excluded_loci} marker STR tidak cocok dan sangat mustahil diturunkan dari pria yang diuji. Secara definitif dan absolut pria ini BUKAN ayah biologis."
    elif prob_paternity > 99.99:
        interp = "PATERNITAS TERBUKTI SECARA MEDIS"
        layman = f"Level kecocokannya mencapai {prob_paternity:.6f}%. Menurut standar pengujian DNA, secara medis pria yang diuji SAH merupakah ayah biologis."
    elif prob_paternity > 99.0:
        interp = "KEMUNGKINAN KUAT PATERNITAS"
        layman = f"Probabilitas berada di angka {prob_paternity:.2f}%. Hampir pasti pria ini ayah kandung, namun dianjurkan untuk memeriksa marker tambahan agar persentasenya melampaui 99.99% untuk kebutuhan legalitas."
    elif excluded_loci > 0 and excluded_loci < 3:
        interp = "INKONKLUSIF (MUTASI TUNGGAL/GANDA)"
        layman = f"Ditemukan {excluded_loci} ketidakcocokan. Hal ini bisa terjadi karena mutasi langka pada sperma ayah, tetapi jumlah data ini terlalu ambigu. Membutuhkan pengujian panel diperluas secara medis."
    else:
        interp = "INKONKLUSIF (DATA MINIM)"
        layman = f"Dengan probabilitas hanya {prob_paternity:.2f}%, kita tidak bisa memastikan kebenaran genetik. Kemungkinan marker yang Anda masukkan sangat sedikit (disarankan uji setidaknya 15+ marker)."

    return cpi, prob_paternity, total_loci, matching_loci, excluded_loci, interp, layman, pd.DataFrame(results)
