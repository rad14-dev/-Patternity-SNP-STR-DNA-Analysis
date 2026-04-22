# Paternity & SNP Analysis Web Dashboard

Aplikasi web interaktif ini akan dibangun menggunakan **Streamlit** untuk mempermudah perhitungan dan visuliasiasi dari skrip `uji_paternitas.py` dan `Single Nuclotide Simulation.py`. Pengguna tidak perlu lagi menggunakan command line; mereka cukup mengunggah file CSV SNP dan langsung melihat hasilnya dalam bentuk grafik dan tabel yang menarik.

## Fitur Utama
1. **Upload Data Interaktif**: Mendukung unggahan file CSV berformat `SNP_ID, Kromosom, Posisi, Mother, Alleged_Father, Child`.
2. **Validasi Kualitas Data (HWE Check)**: Aplikasi akan melakukan penyaringan otomatis pada SNP yang tidak memenuhi ekuilibrium Hardy-Weinberg.
3. **Kalkulator Likelihood Paternitas**: Menghitung probabilitas (Likelihood Ratio) untuk menyimpulkan status paternitas.
4. **Dashboard Visual**: 
   - Tabel interaktif dari data SNP.
   - *Pie chart* perbandingan SNP yang valid dan tidak valid.
   - *Bar chart* atau indikator visual (*gauge chart*) untuk menampilkan rasio Likelihood dengan jelas ke pengguna.

## User Review Required

> [!IMPORTANT]  
> Apakah Anda ingin menggunakan **Streamlit** sebagai framework aplikasinya? 
> Streamlit adalah pilihan terbaik dan paling populer saat ini untuk mengubah skrip data science/bioinformatika Python menjadi aplikasi web interaktif dengan cepat.

## Proposed Changes

### Paternity-SNP-WebApp (Direktori Proyek Baru)
Kita akan membuat folder baru `f:\Python Webdev Workspace\Python for Research\Paternity-SNP-WebApp` untuk memisahkan proyek aplikasi web ini.

#### [NEW] `app.py`
Skrip utama Streamlit yang merender antarmuka pengguna (UI), elemen unggahan file (`st.file_uploader`), serta menampilkan metrik dan grafik menggunakan `plotly` atau `st.bar_chart`.

#### [NEW] `paternity_logic.py`
Ini adalah hasil *refactor* (pembersihan dan penyesuaian) dari skrip `uji_paternitas.py` yang Anda miliki, mengubah fungsi-fungsi perhitungannya agar dapat di-_import_ dengan mudah ke dalam `app.py`.

#### [NEW] `requirements.txt`
Daftar dependensi yang dibutuhkan (seperti `streamlit`, `pandas`, `plotly`) untuk menjalankan aplikasi ini.

## Open Questions

1. Apakah Anda ingin menambahkan kemampuan ekspor (misalnya, pengguna bisa mengunduh laporan *PDF/CSV* dari hasil pemeriksaannya)?
2. Apakah format CSV input (`SNP_ID, Kromosom, Posisi, Mother, Alleged_Father, Child`) akan selalu tetap seperti ini?

## Verification Plan

### Manual Verification
1. Menginstal dependensi melalui command line: `pip install -r requirements.txt`.
2. Menjalankan aplikasi secara lokal: `streamlit run app.py`.
3. Membuka URL lokal di browser Anda.
4. Mengunggah file `snp_data.csv` yang Anda punya, lalu memastikan bahwa metrik Likelihood Ratio yang ditampilkan cocok dengan hasil yang sebelumnya dikeluarkan oleh command line.
