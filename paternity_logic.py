import math
import pandas as pd
from typing import Dict, Tuple, List

# Konstanta
ALLELE_PROBABILITY = 0.5
MOTHER_ALLELE_CONTRIBUTION = 0.5
HWE_PVALUE_THRESHOLD = 0.05

def uji_hwe(genotip_populasi: List[str]) -> Tuple[float, bool]:
    """
    Melakukan uji Hardy-Weinberg Equilibrium menggunakan uji chi-square dengan koreksi Yates.
    Mendukung alel apa pun (tidak hanya A/B) dengan mapping otomatis.
    Return: (p-value, status HWE)
    """
    # Identifikasi alel yang ada dalam populasi
    alel_unik = set()
    for gt in genotip_populasi:
        for char in gt:
            alel_unik.add(char)
    
    alel_list = sorted(list(alel_unik))
    if len(alel_list) > 2:
        # HWE untuk >2 alel lebih kompleks, di sini kita sederhana saja
        return (0.0, False)
    if len(alel_list) < 2:
        return (1.0, True) # Monomorfisme biasanya dianggap HWE
        
    alel1, alel2 = alel_list[0], alel_list[1]
    g11 = alel1 + alel1
    g12 = "".join(sorted(alel1 + alel2))
    g22 = alel2 + alel2
    
    hitungan_genotip = {g11: 0, g12: 0, g22: 0}
    
    for gt in genotip_populasi:
        sorted_gt = "".join(sorted(gt))
        if sorted_gt in hitungan_genotip:
            hitungan_genotip[sorted_gt] += 1
    
    total_individu = len(genotip_populasi)
    if total_individu < 1:
        return (1.0, True)
    
    # Hitung frekuensi alel
    hitung_1 = (2 * hitungan_genotip[g11]) + hitungan_genotip[g12]
    hitung_2 = (2 * hitungan_genotip[g22]) + hitungan_genotip[g12]
    total_alel = hitung_1 + hitung_2
    
    if total_alel == 0:
        return (1.0, True)
    
    p = hitung_1 / total_alel
    q = 1 - p
    
    # Hitung ekspektasi HWE
    eksp_11 = (p ** 2) * total_individu
    eksp_12 = (2 * p * q) * total_individu
    eksp_22 = (q ** 2) * total_individu
    
    # Chi-square calculation dengan koreksi Yates
    chi_square = 0.0
    for obs, exp in zip(
        [hitungan_genotip[g11], hitungan_genotip[g12], hitungan_genotip[g22]],
        [eksp_11, eksp_12, eksp_22]
    ):
        if exp > 0:
            chi_square += ((abs(obs - exp) - 0.5) ** 2) / exp
            
    # Hitung p-value (derajat kebebasan = 1)
    p_value = math.exp(-chi_square / 2)
    status_hwe = p_value > HWE_PVALUE_THRESHOLD
    
    return (p_value, status_hwe)

def parse_vcf(file_path: str) -> pd.DataFrame:
    """
    Membaca file VCF dan mengekstrak genotipe untuk CHILD, MOTHER, FATHER.
    Mengonversi 0/0, 0/1, 1/1 menjadi representasi basa (misal AA, AG, GG).
    """
    rows = []
    header = None
    with open(file_path, "r") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.strip().split("\t")
                continue
            if header:
                data = line.strip().split("\t")
                row_dict = dict(zip(header, data))
                
                # Cari kolom target
                samples = ["CHILD", "MOTHER", "FATHER"]
                result_row = {"SNP_ID": row_dict.get("ID", f"{row_dict['#CHROM']}:{row_dict['POS']}")}
                
                ref = row_dict["REF"]
                alt = row_dict["ALT"].split(",")[0] # Ambil alel alternatif pertama
                
                valid_all = True
                for s in samples:
                    if s in row_dict:
                        gt_info = row_dict[s].split(":")[0] # Ambil GT field
                        if "/" in gt_info:
                            indices = gt_info.split("/")
                        elif "|" in gt_info:
                            indices = gt_info.split("|")
                        else:
                            valid_all = False
                            break
                            
                        if "." in indices:
                            valid_all = False
                            break
                            
                        # Konversi indeks ke basa
                        genotype = ""
                        for idx in indices:
                            if idx == "0": genotype += ref
                            elif idx == "1": genotype += alt
                            else: genotype += "?" # Multi-allelic NOT supported yet
                        result_row[s.capitalize()] = genotype
                    else:
                        valid_all = False
                
                if valid_all:
                    # Rename Father to Alleged_Father to match internal logic
                    result_row["Alleged_Father"] = result_row.pop("Father")
                    rows.append(result_row)
                    
    return pd.DataFrame(rows)

def calculate_allele_probability(
    child_allele: str, father_allele1: str, father_allele2: str
) -> float:
    """Menghitung probabilitas kontribusi alel ayah"""
    if child_allele == father_allele1 or child_allele == father_allele2:
        return 0.5  # Probabilitas 50% untuk setiap alel ayah
    return 0.0

def calculate_snp_likelihood(
    child_genotype: str, mother_genotype: str, father_genotype: str
) -> Tuple[float, float]:
    """Menghitung likelihood untuk satu SNP"""
    # Normalisasi genotipe
    child = ''.join(sorted(child_genotype))
    mother = ''.join(sorted(mother_genotype))
    father = ''.join(sorted(father_genotype))
    
    if len(child) != 2 or len(mother) != 2 or len(father) != 2:
        return 1.0, 1.0
    
    # Hitung probabilitas alel dari orang tua
    prob_paternitas = 0.0
    prob_non_paternitas = 0.0
    
    # Kemungkinan kombinasi alel ibu
    for m_allele in [mother[0], mother[1]]:
        # Alel yang harus diberikan ayah
        required_alleles = [a for a in child if a != m_allele]
        if len(required_alleles) == 0:
            # Kasus dimana alel ibu sudah mencukupi
            if required_alleles:
                prob_paternitas += 0.5 * calculate_allele_probability(required_alleles[0], father[0], father[1])
            else:
                prob_paternitas += 0.5 * 1.0
        else:
            for req_allele in required_alleles:
                prob_paternitas += 0.5 * calculate_allele_probability(req_allele, father[0], father[1])
        
        # Probabilitas non-paternitas (acak populasi)
        prob_non_paternitas += 0.5 * ALLELE_PROBABILITY
    
    return prob_paternitas, prob_non_paternitas

def hitung_likelihood(
    genotipe_anak: Dict[str, str],
    genotipe_ibu: Dict[str, str],
    genotipe_ayah_dugaan: Dict[str, str],
) -> Tuple[float, float]:
    """
    Menghitung likelihood paternitas dengan mempertimbangkan HWE
    """
    likelihood_paternitas = 1.0
    likelihood_bukan_ayah = 1.0

    for snp_id in genotipe_anak:
        paternity_snp, non_paternity_snp = calculate_snp_likelihood(
            genotipe_anak[snp_id],
            genotipe_ibu[snp_id],
            genotipe_ayah_dugaan[snp_id]
        )
        
        likelihood_paternitas *= paternity_snp
        likelihood_bukan_ayah *= non_paternity_snp

    return likelihood_paternitas, likelihood_bukan_ayah

def interpretasi_likelihood_ratio(likelihood_ratio: float) -> Tuple[str, str]:
    """Interpretasi hasil likelihood ratio"""
    if likelihood_ratio > 1000:
        return ("Bukti paternitas SANGAT KUAT", "Kecocokan genetik antara pria yang diuji dan anak sangat tinggi. Terdapat keyakinan yang sangat kuat bahwa pria tersebut adalah ayah biologis yang sah dari anak.")
    elif likelihood_ratio > 100:
        return ("Bukti paternitas KUAT", "Secara statistik, pria yang diuji sangat mungkin adalah ayah biologis anak tersebut, karena probabilitas kecocokan genetiknya kuat.")
    elif likelihood_ratio > 10:
        return ("Bukti paternitas CUKUP", "Ada kecocokan genetik yang mendukung bahwa pria yang diuji adalah ayah dari si anak, namun sangat disarankan untuk menguji lebih banyak marker (SNP) guna mendapatkan hasil yang lebih meyakinkan.")
    elif likelihood_ratio > 1:
        return ("Bukti paternitas LEMAH", "Data genetik secara lemah mengisyaratkan kemungkinan ayah biologis, namun sampel saat ini sama sekali belum cukup untuk membuat kesimpulan yang akurat.")
    elif likelihood_ratio == 1:
        return ("Tidak ada bukti yang cukup", "Berdasarkan data yang diuji, probabilitas menjadi ayah sama persis dengan probabilitas bukan ayah. Kita tidak dapat menarik kesimpulan apa pun. Diperlukan pengujian tambahan yang lebih luas.")
    else:
        return ("Bukti mendukung NON-paternitas", "Data genetik anak memiliki ketidakcocokan yang tidak mungkin diturunkan dari pria yang diuji. Oleh karena itu, pria tersebut terbukti BUKAN ayah biologis dari anak tersebut.")

def process_dataframe(df: pd.DataFrame) -> Tuple[Dict, Dict, Dict, List[str], float, float, float, str, str]:
    """
    Processing SNP data from a pandas DataFrame.
    """
    genotipe_ibu = {}
    genotipe_ayah_dugaan = {}
    genotipe_anak = {}
    snps_tidak_hwe = []

    for _, row in df.iterrows():
        try:
            snp_id = str(row["SNP_ID"]).strip()
        except KeyError:
            continue
            
        genotipe_ibu[snp_id] = str(row.get("Mother", "")).strip().upper()
        genotipe_ayah_dugaan[snp_id] = str(row.get("Alleged_Father", "")).strip().upper()
        genotipe_anak[snp_id] = str(row.get("Child", "")).strip().upper()

        # Periksa HWE untuk SNP ini
        populasi = [genotipe_ibu[snp_id], genotipe_ayah_dugaan[snp_id]]
        _, status_hwe = uji_hwe(populasi)
        if not status_hwe:
            snps_tidak_hwe.append(snp_id)
            
    snp_valid = [snp for snp in genotipe_anak if snp not in snps_tidak_hwe]
    
    ibu_hwe = {snp: genotipe_ibu[snp] for snp in snp_valid}
    ayah_hwe = {snp: genotipe_ayah_dugaan[snp] for snp in snp_valid}
    anak_hwe = {snp: genotipe_anak[snp] for snp in snp_valid}
    
    L_paternitas, L_non = hitung_likelihood(anak_hwe, ibu_hwe, ayah_hwe)
    LR = L_paternitas / L_non if L_non != 0 else float('inf')
    interpretasi, penjelasan_awam = interpretasi_likelihood_ratio(LR)
    
    return genotipe_ibu, genotipe_ayah_dugaan, genotipe_anak, snps_tidak_hwe, L_paternitas, L_non, LR, interpretasi, penjelasan_awam
