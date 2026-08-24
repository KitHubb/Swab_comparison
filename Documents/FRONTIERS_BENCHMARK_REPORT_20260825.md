# Frontiers in Microbiomes 투고용 분석 benchmark 및 비판적 평가

- 작성일: 2026-08-25
- 목표 저널: **Frontiers in Microbiomes**
- 분석 기준: decontam prevalence threshold 0.5 + singleton 제거
- Rarefaction depth: 5,876 reads
- 모든 permutation 기반 검정: 9,999 permutations
- Random seed: 42
- 입력 객체: `Phyloseq/phy_F270R240_260812_v3.rds`, `Phyloseq/phy_F270R240_260812_v4.rds`

## 1. 결론부터

현재 연구는 **Frontiers in Microbiomes에 투고를 준비할 수 있는 수준의 분석적 근거를 확보했다.** 저널의 현재 Impact Factor는 3.0이고 CiteScore는 3.2이며, Methods와 Original Research를 모두 받는다. 다만 이 저널은 단순 관찰적 microbiome survey보다 질문 또는 가설 중심 연구를 요구한다. 따라서 원고의 중심은 “다섯 swab에서 유의한 차이가 없었다”가 아니라 다음 질문이어야 한다.

> 저생체량 피부 16S 자료에서 참여자와 해부학적 부위의 신호가 다섯 가지 integrated swab–medium sampling system과 합리적인 decontamination 및 rare-feature filtering 조건을 바꾸어도 유지되는가?

이번 추가 분석은 이 질문에 대체로 긍정적인 근거를 제공한다. 표본 손실이 없는 decontam threshold 0.1–0.5 범위에서 participant 효과가 가장 컸고 sampling-system 효과는 작고 비유의였다. Leave-one-participant-out 분석에서도 특정 참여자 한 명이 이 결론을 만들지 않았다. 같은 participant×site block 안에서 sampling system을 달리한 거리도 participant나 site가 달라진 비교보다 작았다.

게재 가능성을 제한하는 문제도 남아 있다. 참여자가 5명뿐이며, sampling 위치·순서의 무작위화 여부가 확인되지 않았고, air control 이외의 extraction blank, PCR no-template control, unused-swab control 및 positive control 수행 여부가 불명확하다. qPCR 또는 spike-in 자료가 없어 절대 bacterial biomass를 비교할 수 없다. 이 문제들은 추가 통계로 제거할 수 없으며 Methods와 Limitations에 정확히 기록해야 한다.

## 2. 가장 중요한 기술적 감사 결과: UniFrac tree 교정

기존 post-filter phylogenetic tree는 rooted 상태이고 ASV tip과 정확히 일치했지만, feature pruning 후 5개의 polytomy가 남아 strictly bifurcating tree가 아니었다. `phyloseq`의 fast UniFrac 계산은 이진 tree를 가정하기 때문에 다음 경고가 반복되었다.

> data length is not a sub-multiple or multiple of the number of rows

이를 무시하면 Weighted 및 Unweighted UniFrac 거리가 잘못 계산될 수 있다. 이번 분석에서는 원래 branch length와 tip을 변경하지 않고 `ape::multi2di(random = FALSE)`로 5개의 zero-length edge를 추가하였다.

| Tree audit | 교정 전 | 교정 후 |
|---|---:|---:|
| Tips | 1,689 | 1,689 |
| Internal nodes | 1,683 | 1,688 |
| Edges | 3,371 | 3,376 |
| Zero-length edges | 113 | 118 |
| Rooted | Yes | Yes |
| Strictly binary | No | Yes |
| Tip set preserved | - | Yes |

이 교정 후 UniFrac 분석에서는 tree warning이 발생하지 않았다. **기존 원고와 Figure에 기록된 Weighted/Unweighted UniFrac PERMANOVA 수치는 더 이상 최종값으로 사용하면 안 된다.** Figure 1, Figure S2, Figure S3 및 UniFrac 기반 본문·표를 corrected tree 기준으로 다시 생성해야 한다. Bray–Curtis 결과는 tree의 영향을 받지 않아 기존과 동일하다.

이 문제를 발견하고 재현 가능한 교정 절차와 tree audit 표를 남긴 것은 논문의 분석 신뢰도를 높인다. Methods에는 pruning 후 polytomy를 zero-length branches로 결정론적으로 해소했다는 내용을 짧게 명시한다.

## 3. Canonical singleton 결과

Canonical 조건은 50개 true sample, 1,780 ASVs, 847,074 reads를 보존하였다. Singleton 제거로 손실된 것은 75 ASVs와 75 reads이며, post-decontam read의 99.991%가 유지되었다.

Corrected-tree 단변량 PERMANOVA 결과는 다음과 같다.

| Distance | Participant R² | Site R² | Sampling system R² | System p | System PERMDISP p |
|---|---:|---:|---:|---:|---:|
| Weighted UniFrac | 0.4482 | 0.1682 | 0.0152 | 0.9981 | 0.9898 |
| Unweighted UniFrac | 0.2620 | 0.0633 | 0.0558 | 0.8324 | 0.9519 |
| Bray–Curtis | 0.4080 | 0.1155 | 0.0300 | 0.9947 | 0.7993 |

Participant와 site는 세 거리에서 모두 p=0.0001이었다. Sampling system은 모두 비유의였고 multivariate dispersion 차이도 없었다. Weighted UniFrac에서 participant 효과는 system 효과의 약 29.6배, site 효과는 약 11.1배였다. Bray–Curtis에서는 각각 약 13.6배와 3.9배였다.

Unweighted UniFrac에서는 site R²=0.0633과 system R²=0.0558이 가깝다. 따라서 “모든 거리에서 site 효과가 system 효과보다 압도적으로 크다”고 쓰면 과장이다. 더 안전한 결론은 participant 효과가 일관되게 가장 컸으며, sampling system의 유의한 효과는 검출되지 않았다는 것이다.

## 4. Decontam threshold × rare-feature filtering benchmark

다음 27개 조건을 분석하였다.

- Decontam prevalence threshold: 0.1–0.9
- No additional filtering
- Singleton removal
- Total count ≤10 제거

Threshold 0.5/no-filter 재구성은 저장된 v4와 정확히 일치했다.

- Reconstructed: 1,855 ASVs, 847,149 reads, 50 samples
- Saved v4: 1,855 ASVs, 847,149 reads, 50 samples
- ASV identity: exact match

### 4.1 표본 손실이 없는 범위

Threshold 0.1–0.5의 15개 조합에서는 50개 sample이 모두 rarefaction에 포함되었다.

| Distance | Participant R² range | Site R² range | System R² range | Minimum system p | Minimum PERMDISP p |
|---|---:|---:|---:|---:|---:|
| Weighted UniFrac | 0.4478–0.4592 | 0.1491–0.1690 | 0.0106–0.0152 | 0.9979 | 0.9887 |
| Unweighted UniFrac | 0.2530–0.2985 | 0.0591–0.0649 | 0.0556–0.0626 | 0.4258 | 0.3546 |
| Bray–Curtis | 0.4076–0.4124 | 0.1117–0.1160 | 0.0286–0.0300 | 0.9945 | 0.7993 |

이 범위에서는 결과가 매우 안정적이다. 특히 Weighted UniFrac과 Bray–Curtis의 효과크기는 threshold와 singleton/count≤10 선택에 거의 영향을 받지 않았다.

### 4.2 과도한 decontam의 경계조건

Threshold 0.6부터 대규모 sample loss가 발생했다.

| Threshold | True-sample reads retained | Samples below 5,876 | Rarefied samples |
|---|---:|---:|---:|
| 0.5 | 91.64% | 0 | 50 |
| 0.6 | 49.13% | 17 | 33 |
| 0.7 | 43.15% | 18 | 32 |
| 0.8 | 35.37% | 20 | 30 |
| 0.9 | 20.40% | 36 | 14 |

Threshold 0.9에서 system R²가 겉보기에는 최대 0.33까지 상승했지만 14개 sample만 남은 결과이며 system p는 비유의였다. 이 값은 sampling-system 효과 증가로 해석할 수 없다. 분석 대상이 크게 바뀐 결과이다. Threshold 0.9의 일부 Unweighted UniFrac 조건에서는 기존 non-binary 계산과 달리 corrected tree에서 PERMDISP 유의성은 나타나지 않았지만, sample loss 자체로 이미 본 분석 비교 범위를 벗어난다.

논문에서는 threshold 0.1–0.5를 결론 안정성 범위로 제시하고, 0.6–0.9는 과도한 contamination removal이 데이터와 표본을 손실시키는 boundary condition으로 제시하는 것이 좋다.

## 5. Leave-one-participant-out 결과

각 참여자를 한 명씩 제외하여 5회 반복하였다. 모든 PERMANOVA는 9,999 permutations와 seed 42를 사용했다.

| Distance | Variable | R² range | p-value range | Significant runs |
|---|---|---:|---:|---:|
| Weighted UniFrac | Participant | 0.4144–0.4719 | 0.0001–0.0001 | 5/5 |
| Weighted UniFrac | Site | 0.1572–0.2421 | 0.0001–0.0001 | 5/5 |
| Weighted UniFrac | System | 0.0133–0.0262 | 0.9767–0.9997 | 0/5 |
| Unweighted UniFrac | Participant | 0.1569–0.2785 | 0.0001–0.0007 | 5/5 |
| Unweighted UniFrac | Site | 0.0485–0.1212 | 0.0001–0.0051 | 5/5 |
| Unweighted UniFrac | System | 0.0668–0.0764 | 0.8265–0.8671 | 0/5 |
| Bray–Curtis | Participant | 0.3645–0.4319 | 0.0001–0.0001 | 5/5 |
| Bray–Curtis | Site | 0.1174–0.1506 | 0.0001–0.0001 | 5/5 |
| Bray–Curtis | System | 0.0333–0.0408 | 0.9815–0.9969 | 0/5 |

Sampling-system PERMDISP도 모든 leave-one-out 분석에서 비유의였다. p-value 범위는 Weighted UniFrac 0.1144–0.9874, Unweighted UniFrac 0.7270–0.9752, Bray–Curtis 0.6564–0.8711이었다.

이 결과는 특정 참여자 한 명이 system 비유의 결과를 만든 것이 아니라는 점을 보여준다. 작은 n의 약점을 제거하지는 못하지만, influence diagnostic으로 논문 가치가 높다.

## 6. Participant×site matched-block benchmark

Aggregate UpSet은 다른 participant 또는 site에서 검출된 feature도 system 공통으로 계산할 수 있다. 이를 보완하기 위해 5 participants × 2 sites의 10개 block을 만들고 같은 block 안의 다섯 system을 직접 비교하였다.

### 6.1 Beta-diversity 거리

| Distance | Within participant-site, different systems | Same participant, different sites | Different participants, same site/system |
|---|---:|---:|---:|
| Weighted UniFrac median | 0.1892 | 0.5171 | 0.5301 |
| Unweighted UniFrac median | 0.5563 | 0.6900 | 0.7111 |
| Bray–Curtis median | 0.4196 | 0.7287 | 0.8342 |

같은 participant와 site에서 system만 바꾼 거리가 세 지표 모두 가장 작았다. Weighted UniFrac의 within-block median은 between-participant median의 약 36%, Bray–Curtis는 약 50%, Unweighted UniFrac은 약 78%였다. 이는 abundant community structure의 system 간 재현성이 rare-feature presence/absence보다 높다는 뜻이다.

이 pairwise distance들은 서로 독립적이지 않으므로 일반적인 Wilcoxon 검정을 적용하지 않았다. 본 분석은 효과크기와 분포를 보여주는 descriptive matched benchmark로 사용하고, 공식 추론은 participant-restricted PERMANOVA에 둔다.

### 6.2 Block별 공유 ASV와 Genus

| Level | Median shared feature count | Median intersection/union | Median shared-read fraction | Range of shared-read fraction |
|---|---:|---:|---:|---:|
| ASV | 7.5 | 0.054 | 79.67% | 46.22–98.08% |
| Genus | 7 | 0.116 | 95.12% | 81.26–99.78% |

Richness 기준의 완전 교집합은 작지만, 다섯 system 모두에서 검출된 Genus가 각 block read의 중앙값 95.1%를 차지했다. Aggregate Genus 공유 read fraction은 97.36%였고 system-exclusive Genus read는 0.41%였다.

ASV shared-read fraction은 일부 block에서 46.2%까지 내려갔다. 따라서 ASV 수준까지 모든 system이 동일하다고 주장해서는 안 된다. 논문의 주 해석 단위를 Genus로 두고 V1–V3 species annotation을 putative로 제한하는 현재 방향이 타당하다.

## 7. 반복측정 alpha diversity와 read depth

기존 독립표본 Wilcoxon rank-sum 대신 participant를 block으로 둔 Friedman omnibus test를 수행하였다.

| Site | Outcome | Friedman p |
|---|---|---:|
| Antecubital | Shannon | 0.2052 |
| Antecubital | Observed ASVs | 0.1348 |
| Forehead | Shannon | 0.8911 |
| Forehead | Observed ASVs | 0.6339 |
| Antecubital | Post-filter reads | 0.4510 |
| Forehead | Post-filter reads | 0.8372 |

Sampling system의 omnibus 차이는 검출되지 않았다. 이 결과를 본문 통계로 사용하고, paired Wilcoxon signed-rank 사후 비교는 보조표로 둘 수 있다. n=5에서는 비유의 결과가 equivalence를 증명하지 않는다는 문장을 유지해야 한다.

## 8. 기존 연구와 비교한 위치

### 8.1 Bjerre et al., Scientific Reports 2019

9명의 성인, 165 samples, swab과 scrape, 12개 extraction kit, 16S 및 일부 shotgun을 비교했다. eSwab과 scrape에서 99.3%의 sequence overlap을 보고했고 extraction kit에 따른 library success 차이를 제시했다. 본 연구는 참여자 수와 multi-omics 범위에서는 약하지만, 다섯 integrated sampling system, 두 해부학적 부위, decontam/filtering benchmark 및 repeated-measures effect partitioning에서 차별화된다.

- DOI: https://doi.org/10.1038/s41598-019-53599-z

### 8.2 Ogai et al., Frontiers in Microbiology 2018

7명에서 swab과 tape stripping을 sequence 및 culture로 비교했다. 본 연구에는 culture validation이 없지만, 더 많은 swab–medium system과 bioinformatic robustness 분석이 있다.

- DOI: https://doi.org/10.3389/fmicb.2018.02362

### 8.3 Balacco et al., Frontiers in Microbiomes 2025

16명에서 cotton/eSwab, saline/PBS, swabbing duration 및 storage temperature를 비교하였다. 참여자 수와 실험 요인 분해 측면에서 이 연구가 강하다. 본 연구의 swab material과 medium은 완전 요인설계가 아니므로 독립 효과를 주장할 수 없다. 반면 본 연구는 두 anatomical sites, five integrated systems, contamination threshold 및 rare-feature sensitivity, leave-one-out과 block-level concordance를 함께 제시할 수 있다.

- Frontiers article: https://www.frontiersin.org/journals/microbiomes/articles/10.3389/frmbi.2025.1559981/full

### 8.4 최근 low-biomass method 연구

최근 연구들은 collection negative, extraction negative 및 mock community를 함께 사용하고 qPCR 또는 shotgun 결과와 비교한다. 본 연구는 이러한 validation이 없어 “가장 정확한 sampling system”이나 absolute recovery를 평가할 수 없다. 대신 동일 자료에서 contamination/filtering 결론 안정성을 투명하게 보여주는 방향으로 범위를 제한한다.

- Example: https://pmc.ncbi.nlm.nih.gov/articles/PMC12262292/

## 9. Novelty를 실제로 높이는 요소

다음 네 요소를 한 논리로 연결해야 한다.

1. **Pre-analytical comparison:** 동일 participant와 site에서 다섯 integrated swab–medium systems를 비교하였다.
2. **Biological hierarchy:** participant와 site 신호를 sampling-system 신호와 분리하였다.
3. **Low-biomass analytical robustness:** decontam threshold와 rare-feature filtering이 결론에 미치는 영향을 정량화하였다.
4. **Matched concordance:** aggregate overlap을 넘어 participant×site block 안에서 community distance와 shared-read fraction을 평가하였다.

개별 요소는 기존 문헌에 존재하지만, 이 네 요소를 하나의 반복측정 피부 자료에서 함께 보여주는 것이 본 연구의 실질적 novelty이다.

권장 핵심 문장:

> Across five integrated swab–medium systems, participant- and site-associated bacterial community structure was preserved under decontamination and rare-feature filtering conditions that retained the complete sample set. Concordance was strongest for abundance-weighted and genus-level profiles, whereas ASV-level presence/absence showed greater within-block variability.

## 10. 논문의 한계

1. **Participants n=5:** population-level generalization과 equivalence/non-inferiority 결론이 불가능하다.
2. **Sampling 위치·순서:** 무작위화 또는 rotation 여부가 확인되지 않으면 method와 local spatial heterogeneity가 교락될 수 있다.
3. **Control 구성:** air control 외 extraction blank, PCR NTC, unused-swab 및 positive/mock control 여부가 불명확하다.
4. **Absolute biomass 부재:** total DNA concentration은 bacterial load가 아니며 qPCR/spike-in이 없다.
5. **Non-factorial design:** swab material과 collection medium의 독립 효과를 분리할 수 없다.
6. **Single batch and limited population:** lot, operator, storage duration, sequencing batch 재현성을 평가하지 않았다.
7. **V1–V3 resolution:** species assignment는 putative이며 strain-level 해석이 불가능하다.
8. **ASV concordance heterogeneity:** block별 shared ASV reads가 46–98%로 넓어 Genus 결과보다 일관성이 낮다.
9. **Decontam threshold choice:** 0.5는 data-driven compromise이다. Threshold 0.1–0.5 sensitivity를 함께 보고해 선택 편향을 완화한다.
10. **Tree preprocessing:** 기존 분석의 non-binary tree 문제를 수정했으므로 모든 UniFrac Figure와 manuscript 수치를 재동기화해야 한다.

## 11. 현재 투고 준비도 평가

| 평가 항목 | 현재 수준 | 판단 |
|---|---|---|
| 반복측정 통계 | 보완 완료 | 적절함 |
| PERMANOVA/PERMDISP | 3 distances, 9,999 permutations | 적절함 |
| Influence analysis | Leave-one-participant-out 완료 | 강점 |
| Matched-system concordance | Distance 및 ASV/Genus 완료 | 강점 |
| Decontam/filter benchmark | 27 conditions 완료 | 강점 |
| Reproducibility | Script, seed, sessionInfo, sanity checks | 적절함 |
| Sample size | 5 participants | 주요 약점 |
| Experimental controls | 확인 필요 | 주요 약점 |
| Randomization record | 확인 필요 | 잠재적 치명적 약점 |
| Absolute bacterial load | 없음 | 한계 |
| Final manuscript synchronization | 미완료 | 제출 전 필수 |

**비판적 판정:** 분석 깊이는 Frontiers in Microbiomes의 3점대 방법론 논문을 목표로 할 만하다. 수락 가능성을 결정하는 것은 추가적인 통계 개수가 아니라 실험기록의 투명성, corrected UniFrac 결과와 원고의 완전한 동기화, 작은 표본에서 equivalence를 주장하지 않는 절제된 해석이다. Randomization과 control 기록이 확인되고 원고가 robustness 중심으로 재구성되면 투고 가치가 있다. 해당 기록이 없으면 limitation을 명확히 인정해야 하며, “sampling systems are interchangeable”라는 결론은 피해야 한다.

## 12. 제출 전 필수 작업

1. Corrected binary tree 처리로 `05_Final_main_supplement_figures.Rmd`의 UniFrac Figure와 Table을 재생성한다.
2. 원고의 기존 Weighted UniFrac R²=0.43699 및 site R²=0.26936 등을 corrected 값으로 교체한다.
3. Figure 1, Figure S2, Figure S3 legend에 distance, rarefaction, singleton 및 permutation scheme을 명시한다.
4. Sampling 위치와 채취 순서의 randomization/rotation 기록을 의학·채취 담당자에게 확인한다.
5. Negative, extraction, PCR 및 positive control의 실제 수행 여부를 확인한다.
6. Methods에서 swab material과 medium이 독립 요인이 아님을 명시한다.
7. Main 또는 Supplementary Figure로 threshold 0.1–0.5 effect-size stability와 matched-block distance를 추가한다.
8. Supplementary Table에 leave-one-out, PERMDISP, block agreement 및 tree audit를 포함한다.
9. Abstract와 Discussion의 equivalence/interchangeability 표현을 제거한다.
10. 모든 숫자가 canonical singleton object에서 생성됐는지 최종 audit한다.

## 13. 생성 파일

### Decontam/filter benchmark

- Code: `08_Decontam_filter_benchmark.Rmd`
- Report: `Decontam_filter_benchmark_output/Decontam_filter_benchmark.html`
- Figures: `Decontam_filter_benchmark_output/Figures/`
- Tables: `Decontam_filter_benchmark_output/Tables/`

### Frontiers robustness validation

- Code: `09_Frontiers_robustness_validation.R`
- Figures: `Frontiers_robustness_output/Figures/`
- Tables: `Frontiers_robustness_output/Tables/`
- Reproducibility: `Frontiers_robustness_output/sessionInfo.txt`
- Run metadata: `Frontiers_robustness_output/README_results.txt`

## 14. 저널 정보

Frontiers in Microbiomes는 현재 3.0 Impact Factor, 3.2 CiteScore이며 Web of Science ESCI와 Scopus 등에 색인되어 있다. Methods 및 Original Research article type을 받는다. 단순 observational survey는 범위 밖이라고 명시하므로 본 원고는 hypothesis-driven robustness benchmark로 제출해야 한다.

- Journal scope and metrics: https://www.frontiersin.org/journals/microbiomes/about
- Frontiers metrics: https://www.frontiersin.org/about/impact
