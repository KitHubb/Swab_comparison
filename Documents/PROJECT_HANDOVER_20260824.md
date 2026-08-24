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

## 15. 목표 저널 및 투고 전략 업데이트(2026-08-25)

### 15.1 목표 저널

목표 저널은 **Frontiers in Microbiomes**로 정한다. 현재 자료는 Frontiers in Microbiology에 바로 투고하기에는 표본 수와 대조군 구성, sampling-position 무작위화 여부에서 위험이 크다. Frontiers in Microbiomes에는 피부 저생체량 시료의 채취법 비교와 오염·희귀 feature 처리에 따른 결론의 안정성을 다루는 방법론 연구로 제출하는 방향이 더 적합하다.

논문의 중심 질문은 단순히 “5개 sampling system 사이에 유의한 차이가 없었다”가 아니다. 다음과 같이 설정한다.

> **피부 저생체량 16S 자료에서 참여자와 해부학적 부위의 신호가 다섯 가지 swab–medium sampling system 및 합리적인 오염·희귀 feature 처리 조건을 바꾸어도 유지되는가?**

권장 영문 framing은 다음과 같다.

> Participant- and site-associated community structure remained stable across five integrated swab–medium sampling systems and across reasonable contamination and rare-feature filtering decisions.

이 문장은 sampling system의 동등성이나 완전한 상호교환 가능성을 주장하지 않는다. 표본 수가 5명이므로 “no significant difference”를 “equivalent”로 표현해서는 안 된다.

### 15.2 잠정 제목

우선 권장 제목:

> **Robustness of skin bacterial community profiles across five swab–medium sampling systems and bioinformatic filtering strategies**

대안 제목:

> **Within-subject benchmarking of five swab–medium systems for low-biomass skin microbiome profiling across two anatomical sites**

첫 번째 제목은 decontam 및 filtering sensitivity를 논문의 차별점으로 전면에 배치할 때 사용한다. 두 번째 제목은 실험적 sampling-system 비교를 더 강조할 때 사용한다.

### 15.3 유사 연구와 비교한 차별점

본 연구에서 논문 가치가 있는 부분은 다음과 같다.

1. 동일 참여자에서 5개 integrated swab–medium system을 직접 비교하였다.
2. 미생물 생태가 다른 forehead와 antecubital fossa를 함께 분석하였다.
3. participant, anatomical site 및 sampling system이 설명하는 beta-diversity 변이를 구분하였다.
4. 저생체량 피부 자료에 prevalence 기반 decontam을 적용하고 threshold 0.1–0.9를 검토하였다.
5. singleton, doubleton, count ≤10, prevalence 및 mean-abundance filtering에 따른 데이터 보존과 결과 안정성을 비교하였다.
6. 다섯 sampling system에서 공통 검출된 Genus가 전체 read에서 차지하는 비율을 산출하였다.

유사 swab-method 논문도 소규모 반복측정 설계를 사용한 사례가 있으므로 참여자 5명이라는 이유만으로 연구 가치가 사라지지는 않는다. 다만 본 연구가 기존 문헌보다 강하게 주장할 수 있는 부분은 “어떤 swab이 더 우수하다”가 아니라 **저생체량 시료의 분석 결정이 sampling-system 결론에 얼마나 영향을 주는지 함께 정량화했다는 점**이다.

### 15.4 비판적 한계 및 해석 경계

다음 한계는 제출 전에 원자료 또는 연구기록을 확인하고 원고에 명시해야 한다.

1. **Sampling 위치의 무작위화 또는 순환 배정 여부**
   각 anatomical site에서 Type 1–5의 접촉 위치와 채취 순서가 무작위 또는 참여자별로 순환 배정되었는지 확인해야 한다. 고정된 위치·순서를 사용했다면 sampling system 효과와 피부의 국소 공간 변이가 혼재될 수 있다. 이 문제는 사후 통계로 완전히 보정할 수 없다.

2. **저생체량 연구의 대조군 구성**
   현재 air control을 이용한 decontam 결과가 존재한다. Extraction blank, PCR no-template control, unused-swab control, positive/mock-community control의 수행 여부를 확인해야 한다. 수행하지 않았다면 이를 숨기지 말고 제한점으로 기록한다.

3. **작은 참여자 수**
   참여자는 5명이며 한 참여자 안의 반복 시료가 생물학적으로 독립된 참여자 수를 늘려주지 않는다. 검정력 부족으로 sampling-system 차이를 놓칠 수 있으므로 equivalence 또는 non-inferiority를 주장하지 않는다.

4. **Sampling system, swab material 및 collection medium의 구조적 교락**
   본 설계는 완전 요인설계가 아니다. Swab-tip material과 collection medium의 PERMANOVA는 탐색적 결과로 제시하고 독립적인 인과 효과로 해석하지 않는다.

5. **절대 세균량 자료 부재**
   DNA concentration은 human DNA를 포함한 total DNA이며 bacterial biomass 지표가 아니다. qPCR 또는 spike-in 기반 절대량이 없으므로 relative abundance와 read count 중심의 비교라는 한계를 명시한다.

6. **Aggregate UpSet의 해석 한계**
   전체 자료를 sampling system별로 합친 UpSet은 한 system에서 검출된 feature가 다른 참여자 또는 다른 부위에서 나타나도 공통으로 계산될 수 있다. 따라서 participant×site matched block 내 일치도를 별도로 산출해야 한다.

7. **Species-level 해상도**
   V1–V3 및 SILVA 기반 species assignment는 putative annotation이다. 주 결과는 Genus level로 해석하고 species 결과는 제한적으로 제시한다.

### 15.5 제출 전 필수 추가 분석

#### A. Canonical singleton 분석 확정

- `05_Final_main_supplement_figures.Rmd`에 singleton 제거를 통합한다.
- 분석 시작 시 `taxa_sums(phy) > 1`을 적용하되, post-decontam true-sample 객체에서 적용되는지 확인한다.
- 전체 read는 847,149가 아니라 **847,074 reads**, Genus는 347이 아니라 **344 genera**로 동기화한다.
- Figure 1–3, Figure S1–S3, PERMANOVA, UpSet summary 및 Supplementary Table을 같은 객체에서 재생성한다.

#### B. 반복측정에 맞는 alpha diversity 및 read-count 분석

- 기존 site별 pairwise Wilcoxon rank-sum test는 독립표본 검정이므로 반복측정 구조에 가장 적합하지 않다.
- 우선 선택은 participant를 random intercept로 둔 mixed-effects model이다.
- 표본 수가 작아 모델이 불안정하면 site별 Friedman test를 전체 검정으로 사용하고, 사후 비교에는 paired Wilcoxon signed-rank test와 multiplicity correction을 적용한다.
- DNA concentration은 bacterial biomass 결과로 해석하지 않는다. Post-filter read count는 필요 시 동일 반복측정 구조로 분석한다.

#### C. Beta-diversity 가정과 안정성 검증

- Weighted UniFrac, Unweighted UniFrac 및 Bray–Curtis에 대해 sampling system별 PERMDISP를 수행한다.
- PERMANOVA의 비유의 결과가 dispersion 차이와 혼재되지 않았는지 확인한다.
- Leave-one-participant-out 분석을 수행하여 각 참여자를 한 명씩 제외했을 때 participant, site 및 system의 R²와 p-value 방향이 유지되는지 평가한다.
- 본문에는 단일 p-value보다 R²의 범위와 방향의 안정성을 우선 보고한다.

#### D. Participant×site matched-block 일치도

- 5 participants × 2 sites의 총 10개 matched block을 정의한다.
- 각 block에서 다섯 sampling system 사이의 Weighted UniFrac, Unweighted UniFrac 및 Bray–Curtis 평균 또는 중앙거리(10개 pair)를 산출한다.
- 비교 기준으로 같은 site의 participant 간 거리와 같은 participant의 site 간 거리를 산출한다.
- 가능한 요약 지표:
  - within-block five-system median distance
  - between-participant, same-site median distance
  - within-participant, between-site distance
  - within-block/shared Genus read fraction
  - block별 검출 Genus Jaccard similarity
- 목표는 sampling-system 내 거리보다 participant 또는 site가 달라질 때 거리가 더 커지는지를 직관적으로 제시하는 것이다.

#### E. Decontam threshold × rare-feature filtering 결론 안정성

- Decontam threshold 0.1–0.9 각각에서 downstream 분석 객체를 재구성한다.
- 각 threshold에서 최소한 no additional filter, singleton removal 및 count ≤10 removal을 비교한다.
- 추가 sensitivity 조건으로 prevalence ≥2, prevalence ≥5/10%, mean relative abundance ≥0.1% 및 ≥1%를 유지할 수 있으나, sample loss를 반드시 함께 보고한다.
- 조합별로 다음 값을 한 표에 집계한다.
  - retained ASVs, reads 및 samples
  - rarefaction depth 5,876 미만 sample 수
  - participant, site 및 system의 PERMANOVA R²와 p-value
  - PERMDISP p-value
  - five-system shared ASV/Genus 수와 read fraction
  - system-exclusive read fraction
  - baseline 대비 sample-level composition similarity
- 분석의 목표는 유의한 조건만 선택하는 것이 아니다. 합리적인 처리 범위에서 **participant/site 신호가 system 신호보다 크다는 결론이 유지되는지** 확인하는 것이다.

#### F. Figure 및 Table 권장 구성

- Main Figure 1: sampling system별 alpha/beta diversity. 반복측정 통계와 PERMANOVA strata를 legend에 명시한다.
- Main Figure 2: site별 Genus composition. Mean relative abundance ≥1%인 taxa만 legend에 표시하고 나머지는 `Other`로 합친다.
- Main Figure 3: Genus-level UpSet 또는 matched-block 일치도 Figure. Aggregate UpSet만으로 swab 간 일치성을 주장하지 않는다.
- 신규 Main 또는 Supplementary Figure: decontam/filtering 조건에 따른 retained reads, system R², site R² 및 participant R²의 안정성 plot.
- 신규 Supplementary Figure: leave-one-participant-out PERMANOVA와 PERMDISP 결과.
- 신규 Supplementary Table: 모든 decontam/filtering 조합의 데이터 보존 및 효과크기.

### 15.6 원고에서 즉시 교정할 사항

1. Methods의 “DADA2 without fixed-length truncation (270/240)” 모순을 제거하고 실제 적용한 truncation length를 정확히 기술한다.
2. Singleton 채택 후 total reads, ASV 수, Genus 수, shared-feature 비율을 전부 같은 분석 객체에서 다시 가져온다.
3. 원고의 Figure S4–S7 표기를 실제 최종 출력인 Figure S1–S3과 맞춘다.
4. Figure 1 legend에 Shannon index뿐 아니라 Observed ASVs가 포함됨을 명시한다.
5. Supplementary Table 번호 중복 또는 S2/S3 충돌을 정리한다.
6. Cutadapt version, QIIME 2 version, SILVA release, classifier training 범위 및 confidence 0.7을 실행 기록과 대조한다.
7. Decontam threshold 0.5가 control read 제거와 true-sample 보존의 균형으로 선택되었음을 기술하되, 사후 선택에 따른 편향 가능성을 sensitivity analysis로 보완한다.
8. “sampling systems were equivalent/interchangeable” 대신 “no statistically significant system-associated difference was detected” 또는 “the observed system effect was small relative to participant and site effects”를 사용한다.

### 15.7 투고 판단 기준

다음 조건이 충족되면 Frontiers in Microbiomes 투고본으로 정리한다.

- singleton 기반 canonical 결과와 원고 수치가 완전히 일치한다.
- sampling 위치와 채취 순서의 무작위화 여부가 확인되어 Methods 또는 Limitations에 반영된다.
- 사용한 negative/positive control 종류가 명확히 기록된다.
- 반복측정 통계, PERMDISP 및 leave-one-participant-out 결과가 추가된다.
- matched participant×site block 분석으로 sampling-system 일치도를 직접 제시한다.
- decontam/filtering 변화에도 핵심 효과크기 순서가 유지됨을 보여준다.
- 표본 수가 작다는 점과 equivalence를 증명하지 못한다는 점을 Discussion에 명시한다.

위 조건을 충족하면 본 연구는 단순한 장비·swab 비교를 넘어, **저생체량 피부 microbiome에서 채취 시스템과 bioinformatic filtering이 결과 해석에 미치는 영향을 함께 평가한 재현성 연구**로 제출할 수 있다.

## 16. 업데이트된 실행 순서

1. Sampling 위치·순서 및 실험 대조군 기록 확인
2. Canonical script에 singleton 제거 통합
3. Figure/Table 전체 재생성 및 원고 수치 동기화
4. 반복측정 alpha/read-count 통계 수행
5. PERMDISP 및 leave-one-participant-out PERMANOVA 수행
6. Participant×site matched-block 거리 및 공유 Genus 분석
7. Decontam threshold × filtering sensitivity 확장 분석
8. Frontiers in Microbiomes 형식에 맞춰 Methods, Results, Discussion 및 Supplementary material 정리
9. 최종 원고에서 equivalence 표현, Figure/Table 번호 및 software version 감사
