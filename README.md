# 🌐 Patternity: Web Dashboard Analysis

Bagian ini adalah versi web dari platform Patternity, yang dirancang menggunakan **Streamlit** untuk analisis DNA (SNP & STR) secara interaktif melalui browser.

## 🚀 Cara Menjalankan

### 1. Prasyarat
Pastikan Anda memiliki Python 3.8+ terinstal.

### 2. Persiapan Lingkungan
Disarankan menggunakan virtual environment:
```bash
# Buat venv
python -m venv venv

# Aktifkan venv (Windows)
venv\Scripts\activate

# Aktifkan venv (Linux/Mac)
source venv/bin/activate
```

### 3. Instalasi Dependensi
Instal library yang diperlukan:
```bash
pip install streamlit pandas plotly
```
*(Atau gunakan `pip install -r requirements.txt` jika file tersedia)*

### 4. Jalankan Aplikasi
Jalankan perintah berikut di dalam folder `Web_Version`:
```bash
streamlit run app.py
```

## 🧪 Fitur Analisis
- **STR Analysis**: Unggah file CSV atau GeneMapper (.txt) untuk kalkulasi CPI dan Probability of Paternity.
- **SNP Analysis**: Mendukung format VCF dan CSV dengan uji Hardy-Weinberg Equilibrium (HWE).
- **Interactive Charts**: Visualisasi grafik distribusi genetik dan gauge likelihood ratio.

---
**Catatan Keamanan**: Versi ini menggunakan pemrosesan file sementara yang aman (`tempfile`) dan tidak menyimpan data genetik Anda secara permanen di server.
