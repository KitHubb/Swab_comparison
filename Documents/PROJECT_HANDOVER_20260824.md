# DT Swab 연구 분석 인수인계 보고서

- 작성일: 2026-08-24
- 연구 저장소: [KitHubb/Swab_comparison](https://github.com/KitHubb/Swab_comparison)
- 기본 브랜치: `main`
- 본 보고서 작성 직전 분석 커밋: `c35772e` (`Add microbiome filtering sensitivity analysis`)
- 주 분석 프로젝트: `Rproj_DT_Swab`

## 1. 인수인계 핵심 요약

본 연구는 5가지 피부 swab–collection medium 조합이 피부 세균군집 분석 결과에 미치는 영향을 비교하는 V1–V3 16S rRNA 기반 방법론 연구이다. 분석의 기본 원칙은 **기존 논문의 분석 구조를 유지하면서 prevalence 기반 decontam을 적용한 뒤 전체 통계와 Figure/Table을 재생성하는 것**이다.

현재까지 다음 작업이 완료되었다.

1. HN00182797 데이터셋을 주 분석 데이터로 선정하였다.
2. Cutadapt, QIIME 2 import, DADA2, SILVA taxonomy, phylogenetic tree 생성이 가능한 Nextflow 워크플로를 구축·실행하였다.
3. DADA2 truncation은 forward 270 bp, reverse 240 bp를 사용하였다.
4. `decontam` prevalence 방법의 threshold 0.1–0.9 민감도 분석을 수행하고, threshold 0.5를 본 분석 기준으로 선택하였다.
5. post-decontam phyloseq 객체, metadata, 주요 Figure/Table을 생성하였다.
6. alpha diversity, beta diversity, community composition, PERMANOVA, ASV/Genus 공유 분석을 재수행하였다.
7. singleton, doubleton, count ≤10 제거 결과를 별도로 비교하였다.
8. 추가 filtering sensitivity analysis를 완료하였다.

현재의 최종 분석 권고는 **decontam threshold 0.5 적용 후 singleton(total count = 1 ASV) 제거**이다. Singleton 제거는 전체 read의 99.991%와 50개 true sample을 모두 보존하면서 극저빈도 ASV만 제거하였다.

중요: 이 권고는 민감도 분석에서 확정되었지만, canonical script인 `05_Final_main_supplement_figures.Rmd`에는 아직 singleton 제거가 통합되지 않았다. 다음 작업자는 singleton 조건을 본 분석 script에 반영하고 최종 Figure/Table 및 원고 수치를 다시 동기화해야 한다.

## 2. 연구 설계와 분석 기준

### 2.1 연구 대상

- 대상: 피부 bacterial microbiome
- 증폭 부위: 16S rRNA V1–V3
- 참여자: 5명
- 해부학적 부위: antecubital fossa, forehead
- sampling system: Type 1–Type 5
- true sample: 50개
- negative/control sample: decontam 판정에 포함

### 2.2 Sampling system 매핑

| Sampling system | 원래 변수 `swab_type2` |
|---|---|
| Type 1 | Copan |
| Type 2 | Puritan_gel |
| Type 3 | Puritan_buffer |
| Type 4 | Omni |
| Type 5 | Eswab |

Swab material과 collection medium은 sampling system과 완전히 독립된 factorial design이 아니다. 따라서 이 두 변수의 PERMANOVA 결과는 탐색적 결과로만 해석한다.

### 2.3 현재 확정된 분석 설정

- 주 raw dataset: `HN00182797`
- Cutadapt quality cutoff: Q20
- DADA2 truncation: forward 270 bp, reverse 240 bp
- taxonomy database: SILVA 138.1, 99% classifier
- taxonomy confidence: 0.7
- decontam: prevalence method, threshold 0.5
- rarefaction depth: 5,876 reads
- 추가 rare-feature filtering 권고: singleton 제거
- community composition 주 해석 수준: Genus
- V1–V3 Species annotation: putative assignment로 제한하여 신중하게 해석
- taxonomy 정제: chloroplast, mitochondria, kingdom-level unassigned 및 eukaryote 제외

## 3. Sequence preprocessing 및 Nextflow 작업

### 3.1 서버 경로

- Nextflow pipeline: `/data/software/nextflow/amplicon_16S_v1v3_qiime_nf/`
- 분석 프로젝트: `/data/home2/ksy/260811_DT_swab/`
- 주 raw data: `/data/FASTQ/HN00182797/DT/`
- 비교용 raw data: `/data/FASTQ/HN00194709_2nd_HN00182797/DT/`
- 주 결과: `/data/home2/ksy/260811_DT_swab/Output/HN00182797/`
- work directory: `/data/home2/ksy/260811_DT_swab/work/HN00182797/`
- SILVA classifier: `/data/Reference/QIIME2-2025.7/Bacteria/SILVA/silva-138-99-nb-classifier.qza`
- QIIME 2 Singularity image: `/data/software/singularity/qiime2_amplicon_2025.7`

### 3.2 수행된 분석 흐름

1. Raw FastQC/MultiQC
2. V1–V3 primer 제거 및 Q20 filtering(Cutadapt)
3. Clean FastQC/MultiQC
4. QIIME 2 paired-end import
5. DADA2 denoising, merging, ASV inference 및 chimera removal
6. SILVA taxonomy classification
7. Phylogenetic tree 생성
8. Alpha/beta diversity 산출(필요 시)

### 3.3 주의사항

Windows의 `Nextflow/` 폴더는 현재 `Swab_comparison` Git 저장소에 포함되지 않으며 별도 Git 저장소도 아니다. 다른 컴퓨터로 완전히 이전하려면 이 폴더를 별도로 복사하거나 별도의 pipeline GitHub 저장소에 먼저 올려야 한다.

Raw FASTQ와 서버 work directory는 현재 GitHub에 포함되지 않는다. Nextflow 재실행은 서버 경로와 Singularity image 접근이 가능한 환경에서 수행해야 한다.

## 4. Decontam 분석

### 4.1 수행 내용

`01_Preprocessing.Rmd`에서 `decontam` prevalence 방법을 이용하여 threshold 0.1–0.9를 비교하였다. 기준은 control read 제거와 true-sample 정보 보존의 균형이었다.

### 4.2 Threshold 0.5 결과

| 항목 | 결과 |
|---|---:|
| 제거 contaminant ASV | 107 |
| contaminant relative abundance | 14.39% |
| 전체 read 보존 | 85.610% |
| true-sample read 보존 | 91.644% |
| control read 보존 | 24.322% |
| true-sample ASV 보존 | 94.836% |
| control ASV 보존 | 58.846% |
| Control-status PERMANOVA R² | 0.0572 |
| Control-status PERMANOVA p | 0.0026 |

Threshold 0.5는 true sample의 ASV와 read를 대부분 유지하면서 control read의 약 75.7%를 제거하였다. 이 결과를 현재 본 분석 기준으로 사용한다.

### 4.3 주요 결과 파일

- `Tables/Decontam_summary.txt`
- `Tables/Decontam_permanova_summary.txt`
- `Phyloseq/phy_F270R240_260812_v4.rds`

`phy_F270R240_260812_v4.rds`는 post-decontam 분석의 공통 시작 객체이다.

## 5. Filtering sensitivity analysis

### 5.1 비교한 조건

모든 조건은 post-decontam true-sample phyloseq 객체에 독립적으로 적용하였다.

1. 추가 filtering 없음
2. total count ≤1 제거(singleton 제거)
3. total count ≤2 제거
4. total count ≤10 제거
5. prevalence ≥2 samples
6. prevalence ≥5 samples
7. prevalence ≥10% samples
8. mean relative abundance ≥0.1%
9. mean relative abundance ≥1%

True sample이 50개이므로 prevalence ≥10%는 prevalence ≥5 samples와 정확히 같은 조건이다.

### 5.2 Filtering별 데이터 보존

| Filtering | ASV 보존 | Read 보존 | 5,876 reads 미만 샘플 |
|---|---:|---:|---:|
| 추가 filtering 없음 | 1,855 | 100.000% | 0 |
| Singleton 제거 | 1,780 | 99.991% | 0 |
| Count ≤2 제거 | 1,714 | 99.976% | 0 |
| Count ≤10 제거 | 1,190 | 99.592% | 0 |
| Prevalence ≥2 | 608 | 96.771% | 2 |
| Prevalence ≥5 | 197 | 90.385% | 2 |
| Prevalence ≥10% | 197 | 90.385% | 2 |
| Mean abundance ≥0.1% | 82 | 87.420% | 2 |
| Mean abundance ≥1% | 10 | 67.733% | 13 |

Prevalence와 abundance filtering은 sampling-system R²를 줄일 수 있지만, 동시에 ASV와 sample을 많이 제거한다. 이 경우 sampling-system 효과의 감소가 생물학적 안정성 때문인지 sample loss 때문인지 분리하기 어렵다.

### 5.3 Singleton 제거 결과

- 제거 ASV: 75개
- 제거 read: 75 reads
- 보존 ASV: 1,780개
- 보존 read: 847,074 reads(99.991%)
- 5,876 reads rarefaction에 포함된 sample: 50/50
- unfiltered composition 대비 median similarity: 99.996%
- minimum sample-level similarity: 99.952%

### 5.4 Singleton 기준 단변량 PERMANOVA

| Distance | Variable | R² | F | p-value | Permutation strata |
|---|---|---:|---:|---:|---|
| Weighted UniFrac | Participant identity | 0.43429 | 8.6364 | <0.001 | 없음 |
| Weighted UniFrac | Anatomical site | 0.27424 | 18.1378 | <0.001 | Participant |
| Weighted UniFrac | Sampling system | 0.02444 | 0.2819 | 0.8008 | Participant |
| Unweighted UniFrac | Participant identity | 0.19855 | 2.7871 | <0.001 | 없음 |
| Unweighted UniFrac | Anatomical site | 0.05299 | 2.6858 | <0.001 | Participant |
| Unweighted UniFrac | Sampling system | 0.06153 | 0.7376 | 0.8789 | Participant |
| Bray–Curtis | Participant identity | 0.40797 | 7.7523 | <0.001 | 없음 |
| Bray–Curtis | Anatomical site | 0.11551 | 6.2687 | <0.001 | Participant |
| Bray–Curtis | Sampling system | 0.02996 | 0.3475 | 0.9947 | Participant |

Participant identity는 strata 없이 검정하였다. Anatomical site, sampling system, swab-tip material 및 collection medium은 repeated-measures 구조를 반영하여 participant 내에서 permutation하였다.

### 5.5 Alpha diversity

- Shannon index와 Observed ASVs를 비교하였다.
- Site별 Type 1–Type 5 pairwise Wilcoxon test에 Bonferroni correction을 적용하였다.
- 9개 filtering 조건에서 총 360개 pairwise comparison을 수행하였다.
- Adjusted p<0.05인 sampling-system 비교는 0개였다.

### 5.6 Sampling system 공유 feature

Singleton 제거 후 전체 데이터에서 다음 결과를 얻었다.

| Level | 전체 feature | 5개 Type 공통 feature | 공통 feature가 차지하는 read | System-exclusive read |
|---|---:|---:|---:|---:|
| ASV | 1,780 | 119 | 86.332% | 3.373% |
| Genus | 344 | 59 | 97.404% | 0.342% |

공유 feature의 개수는 전체 richness의 일부이지만, 공통 Genus가 전체 read의 97% 이상을 차지한다. Sampling-system-exclusive feature는 대부분 저풍부도 성분이다. Swab 자체 유래 미생물로 단정하려면 extraction blank, unused-swab control, lot information 및 독립적 반복 검증이 추가로 필요하다.

### 5.7 결론

현재 자료에서는 singleton 제거가 가장 낮은 수준의 추가 filtering이면서 데이터 보존과 재현성을 모두 만족한다. Count ≤2와 count ≤10 제거는 sensitivity analysis로 제시할 수 있다. Prevalence ≥2 이상 또는 mean abundance filtering은 low-biomass sample의 read depth와 sample inclusion을 변경하므로 본 분석의 기본 filtering으로 사용하지 않는 것이 안전하다.

## 6. 현재 Figure 및 Table 작업

### 6.1 완료된 주요 Figure

- `Final_main_supplement_figure_output/Figure1_bacterial_diversity.*`
- `Final_main_supplement_figure_output/Figure2_bacterial_community_composition.*`
- `Final_main_supplement_figure_output/Figure3_UpSet_ASV.*`
- `Final_main_supplement_figure_output/Figure3_UpSet_Genus.*`
- `Final_main_supplement_figure_output/FigureS1_participant_composition.*`
- `Final_main_supplement_figure_output/FigureS2_UniFrac_heatmaps.*`
- `Final_main_supplement_figure_output/FigureS3_participant_alpha_beta.*`

Figure 3 UpSet plot은 ASV 및 Genus level로 생성하였다. Sampling-system label은 45도, intersection bar 위 숫자는 90도로 회전하고 검정색으로 표시하였다. 현재 원고의 주 Figure 3은 Genus-level 결과를 사용하는 방향이다.

### 6.2 Filtering sensitivity 결과

- Code: `06_Filtering_sensitivity_analysis.Rmd`
- HTML: `Filtering_sensitivity_output/Filtering_sensitivity_analysis.html`
- Figures: `Filtering_sensitivity_output/Figures/`
- Tables: `Filtering_sensitivity_output/Tables/`
- 통합 판단표: `Filtering_sensitivity_output/Tables/filtering_decision_summary.csv`

### 6.3 Relative-abundance Table

Figure 2 legend와 Table은 mean relative abundance ≥1%인 Phylum/Genus/Species 목록을 기준으로 재정리하였다. 1% 미만 항목은 `Other`로 합산한다.

주요 파일:

- `Tables/TableS3_Antecubital_relative_abundance_taxa_over_1pct.csv`
- `Tables/TableS4_Forehead_relative_abundance_taxa_over_1pct.csv`
- `Tables/TableS5_Antecubital_pairwise_taxa_over_1pct.csv`
- `Tables/TableS6_Forehead_pairwise_taxa_over_1pct.csv`
- `Tables/Table_S3-S6_taxa_over_1pct.xlsx`

## 7. 원고에서 반영된 핵심 해석

1. Sampling system 간 alpha diversity와 beta diversity의 유의한 차이는 확인되지 않았다.
2. Participant identity와 anatomical site가 sampling system보다 큰 community variation을 설명하였다.
3. Swab material과 collection medium 분석은 study design상 독립 효과로 해석하지 않는다.
4. Community composition은 Genus level을 중심으로 해석한다.
5. V1–V3 Species assignment는 putative annotation으로 제한한다.
6. Low-biomass skin sample 특성상 contamination control과 filtering sensitivity를 명시한다.
7. DNA concentration은 총 DNA 농도이며 bacterial DNA-specific measurement가 아니므로, bacterial biomass의 직접 지표로 해석하지 않는다.
8. 초기 library preparation fail은 재구축 후 sequencing되었으므로, 최종 분석은 성공적으로 sequencing된 library를 사용했다고 기록하고 초기 실패를 sample exclusion처럼 기술하지 않는다.

## 8. GitHub에서 재현 가능한 파일

현재 Git에 포함되어 다른 컴퓨터로 받을 수 있는 주요 입력은 다음과 같다.

- `Phyloseq/phy_F270R240_260812.rds`
- `Phyloseq/phy_F270R240_260812_v3.rds`
- `Phyloseq/phy_F270R240_260812_v4.rds`
- `Phyloseq/phy_F270R240_260812_v4_rarefy_5876.rds`
- `metadata/DT_metadata.tsv`
- Rmd/R scripts
- 최종 Figure/Table
- Filtering sensitivity HTML, Figure 및 CSV
- 문헌 비교 보고서(PDF/DOCX)

`renv.lock`은 사용자의 결정에 따라 포함하지 않았다. 대신 filtering sensitivity output에 `sessionInfo.txt`를 저장하였다.

## 9. GitHub에 포함되지 않는 자료

다음 자료는 현재 `Swab_comparison` 저장소 밖에 있으므로 별도 이전이 필요하다.

### 9.1 원고 및 검토 문서

로컬 경로: `D:\KSY\Project\2.DT_Swab[논문작업pf남경화]\Paper\`

가장 최근 원고:

- `01_CoST_공저자검토용_통합원고_20260726_ksy2.docx`
- `01_CoST_공저자검토용_통합원고_20260726_ksy2.zip`

기타 중요 문서:

- `CoST_manuscript_revision_audit_20260824.docx`
- `04_CoST_재분석_원고_통합점검표_20260824.xlsx`
- `HN00182797_Original_Sample_QC_and_NGS_Library_QC.xlsx`
- `CoST_Supplementary_materials.docx`

### 9.2 Nextflow pipeline

로컬 경로: `D:\KSY\Project\2.DT_Swab[논문작업pf남경화]\Nextflow\`

이 폴더는 현재 Git repository가 아니다. 서버의 pipeline 폴더와 함께 별도 백업이 필요하다.

### 9.3 Raw data 및 서버 결과

FASTQ, Nextflow work directory, QIIME 2 artifact 전체는 GitHub에 포함하지 않았다. 서버 또는 별도 저장장치에서 보존해야 한다.

## 10. 다른 컴퓨터에서 시작하는 방법

### 10.1 신규 clone

```bash
git clone git@github.com:KitHubb/Swab_comparison.git
cd Swab_comparison
git checkout main
git pull origin main
```

SSH가 설정되지 않은 컴퓨터에서는 다음 URL을 사용할 수 있다.

```bash
git clone https://github.com/KitHubb/Swab_comparison.git
```

### 10.2 R 프로젝트 열기

`Rproj_DT_Swab.Rproj`를 RStudio에서 연다. 최근 filtering sensitivity 분석은 다음 환경에서 완료되었다.

- R 4.5.2
- phyloseq 1.54.2
- vegan 2.7-5
- dplyr 1.2.1
- ggplot2 4.0.3
- patchwork 1.3.2
- ComplexUpset
- pheatmap
- RColorBrewer
- ggpubr
- rmarkdown 2.31

전체 버전은 `Filtering_sensitivity_output/sessionInfo.txt`에서 확인한다.

### 10.3 분석 실행 순서

이미 저장된 post-decontam RDS부터 시작할 경우:

1. `06_Filtering_sensitivity_analysis.Rmd` 실행
2. singleton을 본 분석 조건으로 확정
3. `05_Final_main_supplement_figures_singleton.Rmd` 실행
4. singleton 결과를 canonical final output으로 동기화
5. 원고의 Method, Results, Figure legend 및 Supplementary Table 수치 갱신

Decontam부터 완전히 재생성할 경우:

1. `01_Preprocessing.Rmd`
2. `02_Microbiome_analysis.Rmd` 또는 최종 분석에 사용된 수정본
3. `05_Final_main_supplement_figures_singleton.Rmd`
4. `06_Filtering_sensitivity_analysis.Rmd`

`01_Preprocessing.Rmd` 재실행은 phyloseq RDS와 decontam summary를 갱신할 수 있으므로 기존 결과를 보존한 뒤 수행한다.

### 10.4 명령행 RMarkdown 실행 예시

```r
rmarkdown::render(
  "06_Filtering_sensitivity_analysis.Rmd",
  output_file = "Filtering_sensitivity_analysis.html",
  output_dir = "Filtering_sensitivity_output",
  clean = TRUE
)
```

Windows에서 Pandoc을 찾지 못하면 RStudio에서 Knit하거나 `RSTUDIO_PANDOC`을 RStudio의 Pandoc 경로로 지정한다.

## 11. 다음 작업 우선순위

### Priority 1. Singleton을 canonical 본 분석에 반영

- `05_Final_main_supplement_figures.Rmd`에 `taxa_sums(phy) > 1` 조건을 추가하거나 검증된 singleton Rmd를 canonical script로 승격한다.
- Figure 1–3, Figure S1–S3, Table 및 PERMANOVA 결과를 재생성한다.
- 기존 no-filter 결과와 파일명이 섞이지 않도록 output directory를 먼저 구분한 뒤 최종본만 교체한다.

### Priority 2. 원고 수치 동기화

- Methods에 singleton 제거와 rarefaction depth 5,876을 명시한다.
- Results의 PERMANOVA 수치를 singleton 결과로 갱신한다.
- Figure 3의 shared ASV/Genus 비율을 singleton 결과로 갱신한다.
- Filtering sensitivity analysis를 Supplementary Methods/Results 또는 Supplementary Figure로 추가할지 결정한다.

### Priority 3. 최종 원고와 Supplementary material 검수

- 본문 내 Figure/Table 번호와 실제 파일명을 대조한다.
- Figure legend의 filtering level, distance metric, permutation strata를 명시한다.
- Table S3–S6 유지 여부를 최종 결정한다. Pairwise taxon comparison이 주 결론에 기여하지 않으면 축소 또는 제거를 고려한다.
- Taxonomic terminology, `unclassified`, `Other`, `putative species` 표현을 통일한다.

### Priority 4. 외부 파일 이전

- `Paper/` 폴더를 안전한 저장장치 또는 별도 private repository로 이전한다.
- `Nextflow/` 폴더를 별도 Git repository로 만들고 서버 pipeline과 동기화한다.
- Raw FASTQ, QIIME 2 artifacts, classifier 및 server output의 보관 위치를 문서화한다.

## 12. 현재 로컬에만 남아 있는 미커밋 파일

본 보고서 작성 시점에 다음 변경은 기존 작업물로 판단하여 filtering sensitivity 커밋에 포함하지 않았다.

```text
M  05_Final_main_supplement_figures_doubleton.md
M  05_Final_main_supplement_figures_singleton.md
M  05_Final_main_supplement_figures_tenton.md
M  figure/figure-2-1.png
?? Final_main_doubleton/Tables/
?? Final_main_singleton/Tables/
?? Final_main_tenton/Tables/
```

이 파일들은 현재 컴퓨터에만 남아 있을 수 있다. 다른 컴퓨터로 완전히 이전하려면 내용을 검토한 후 별도 커밋하거나 수동 복사해야 한다. 특히 singleton 결과를 최종 분석으로 사용할 예정이므로 `Final_main_singleton/Tables/`는 우선 확인 대상이다.

## 13. 최종 분석 판단

현재 근거는 다음 결론을 지지한다.

- 피부 microbiome variation은 sampling system보다 participant와 anatomical site의 영향을 더 크게 받는다.
- Sampling system 효과는 Weighted UniFrac, Unweighted UniFrac 및 Bray–Curtis에서 유의하지 않다.
- 5개 sampling system이 공유하는 Genus가 전체 community read의 대부분을 차지한다.
- Singleton 제거는 sample과 read를 사실상 모두 유지하면서 극저빈도 ASV의 영향을 줄인다.
- 더 강한 prevalence 또는 abundance filtering은 결과의 겉보기 일치도를 높일 수 있지만 ASV 및 low-depth sample 손실을 동반한다.

따라서 최종 원고의 주 분석은 **decontam threshold 0.5 + singleton 제거 + rarefaction depth 5,876 + SILVA taxonomy + Genus 중심 해석**으로 정리하는 것이 현재 자료와 가장 잘 맞는다.

## 14. 이전 완료 체크리스트

- [ ] 새 컴퓨터에서 `Swab_comparison` clone 또는 pull 완료
- [ ] `Phyloseq/` RDS 4개와 `metadata/DT_metadata.tsv` 확인
- [ ] R package 설치 및 `sessionInfo.txt`와 비교
- [ ] `06_Filtering_sensitivity_analysis.Rmd` 재실행 확인
- [ ] singleton canonical 분석 반영
- [ ] Figure/Table 재생성 및 원고 수치 동기화
- [ ] `Paper/` 폴더 별도 이전
- [ ] `Nextflow/` 폴더 별도 이전 또는 Git 저장소 생성
- [ ] 서버 raw data/output 경로 접근 확인
- [ ] 로컬 미커밋 singleton/doubleton/tenton 결과 처리
