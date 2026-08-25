# Decontam 및 feature-filtering 통합 민감도 분석 보고서

## 분석 목적

본 분석은 다섯 sampling system의 비교라는 논문 주제를 유지하면서, low-biomass 피부 미생물 자료에서 decontam threshold와 희귀 feature 제거 기준이 결론에 미치는 영향을 평가한다. 기존 QIIME 2 계통수는 수정하지 않았으며, 모든 순열 검정은 seed 42와 9,999 permutations를 사용했다.

## 통합 분석 구성

- Decontam 축: 미적용(None)과 prevalence threshold 0.1–0.9를 비교했다. 이 축에는 모든 조건에 동일하게 singleton 제거를 적용했다.
- 후처리 필터 축: decontam threshold 0.5 자료에서 no filter, total count ≤1·2·10 제거, prevalence ≥2·5·10%, mean relative abundance ≥0.1%·1%를 비교했다.
- beta diversity: 기존 tree로 weighted/unweighted UniFrac을 계산하고 Bray–Curtis를 함께 평가했다.
- 반복측정 구조: sampling system 및 anatomical site 검정은 participant 내 제한 순열을 사용했다. Participant identity 검정은 unrestricted permutation을 사용했다.
- 추가 평가: PERMDISP, alpha-diversity Friedman test, sample composition similarity, ASV/genus 공유율, participant×site matched-block 공유율, Procrustes/PROTEST를 포함했다.

## 검증 결과

최종 sanity check 12개가 모두 통과했다. Threshold 0.5 자료는 기존 저장 v4와 ASV 수, true-sample read 수 및 ASV identity가 일치했다. 전처리 phyloseq tree 객체는 분석 전후 동일했고, tree 변환 함수는 사용하지 않았다.

## 주요 결과

Decontam threshold 0.5에서는 true-sample reads의 91.64%가 유지된 반면 negative-control reads는 24.32%만 유지되었다. 50개 true sample이 모두 rarefaction depth 5,876을 충족했다. Threshold 0.1–0.5에서도 모든 50개 sample이 유지됐지만, 0.6에서는 33개, 0.9에서는 14개로 감소했다. 따라서 0.6–0.9 결과는 threshold 자체와 sample loss의 영향을 분리해 해석할 수 없다.

Threshold 0.1–0.5에서 weighted UniFrac sampling-system 효과는 R²=0.0145–0.0244였고 모두 비유의였다. Threshold 0.5 결과는 R²=0.02444, p=0.8027이었다. 같은 조건에서 participant identity는 R²=0.43429, p<0.001, anatomical site는 R²=0.27424, p<0.001이었다. Sampling-system PERMDISP도 세 거리에서 모두 비유의였다(weighted UniFrac p=0.9151; unweighted UniFrac p=0.9418; Bray–Curtis p=0.7993).

Threshold 0.5에서 sampling system 간 Shannon index와 observed ASV 차이는 두 anatomical site 모두 비유의였다. Genus-level에서 다섯 system에 공통인 genera는 read abundance의 97.36%를 차지했다. Decontam 미적용 조건에서도 이 비율은 97.24%로 유사했다.

Singleton, doubleton 및 count ≤10 제거는 read를 각각 99.99%, 99.98%, 99.59% 유지했다. Prevalence 또는 mean-abundance 필터는 더 많은 ASV와 일부 sample을 제거했다. 특히 mean abundance ≥1%는 ASV 10개와 sample 37개만 남겨 sampling-system 비교용 기본 분석으로는 과도하다.

Ordination은 decontam threshold 0.1–0.2에서 미적용 자료와 거의 동일했고, threshold 0.3–0.5에서도 높은 일치도를 유지했다. Weighted UniFrac Procrustes correlation은 threshold 0.5에서 0.9722였다. Threshold 0.6부터 sample loss와 함께 일치도가 감소했다.

## 동일 50개 샘플의 sampling-system 거리 검증

강한 필터에서 sample loss가 발생하는 문제를 피하기 위해, decontam threshold 0.5의 동일한 50개 biological sample을 모든 필터 조건에서 고정하였다. 각 participant×anatomical-site block 안에서 다섯 sampling system 사이의 pairwise distance를 계산하고, participant 간 거리와 anatomical-site 간 거리를 병렬로 평가하였다.

강한 prevalence 및 mean-abundance filtering은 Bray–Curtis와 unweighted UniFrac에서 같은 participant-site 내 sampling-system 거리를 감소시켰다. 예를 들어 prevalence ≥5 samples에서는 sampling-system 거리가 Bray–Curtis에서 19.55%, unweighted UniFrac에서 38.30% 감소했다. Mean relative abundance ≥1%에서는 각각 47.44%와 85.11% 감소했다. Benjamini–Hochberg false discovery rate 보정 후에도 이 감소는 유의했다.

Binary Jaccard와 Aitchison distance에서도 같은 방향이 확인되었다. Prevalence ≥5 samples에서 system 간 거리는 Binary Jaccard 27.28%, Aitchison 37.82% 감소했다. Mean relative abundance ≥1%에서는 각각 65.25%와 81.12% 감소했다. Aitchison은 raw ASV count에 pseudocount 0.5를 더한 뒤 centered log-ratio 변환하여 계산했다. Pseudocount 0.1, 0.5 및 1에서 얻은 거리행렬의 Spearman correlation은 0.9869–0.9993으로 높아, 0.5 선택에 따른 거리 순위 변화는 작았다.

Singleton 및 doubleton 제거에서도 일부 paired test의 BH-adjusted p-value는 0.05 미만이었지만 거리 변화량은 대체로 1% 미만이었다. 따라서 통계적 유의성만으로 filtering의 실질적 효과를 판단하지 않고 median percent change와 biological-signal/system-distance ratio를 함께 해석해야 한다.

이러한 수렴은 모든 거리 지표에서 일관되지 않았다. Weighted UniFrac sampling-system 거리는 prevalence ≥5 samples에서 231.77%, mean relative abundance ≥1%에서 122.45% 증가했다. Participant 및 site 간 weighted UniFrac 거리도 증가했지만, sampling-system 거리의 상대적 증가가 더 커 participant-to-system distance ratio는 no filter 대비 약 41–55% 수준으로 감소했다. 따라서 강한 filtering은 sampling systems의 생물학적 일치도를 보편적으로 개선한 것이 아니라, 거리 지표가 반영하는 community component를 바꾸었다.

Singleton, doubleton 및 total count ≤10 제거에서는 세 거리의 sampling-system 차이가 작았다. Total count ≤10 제거 시 같은 block 내 sampling-system 거리는 Bray–Curtis에서 0.74% 감소하고 unweighted 및 weighted UniFrac에서 각각 1.11%와 5.48% 증가했으며, 보정 후 유의하지 않았다. 이 범위는 희소 feature를 일부 줄이면서 원래 community geometry를 가장 잘 유지한 보수적 filtering 범위로 해석할 수 있다.

## 논문 해석

주 결론은 sampling-system 간 차이가 participant와 anatomical site의 차이보다 작다는 것이다. 이 결론은 decontam 미적용부터 threshold 0.5까지, 그리고 여러 희귀 feature 제거 기준에서 유지됐다. Decontam 분석은 sampling-system 비교를 대체하는 별도 주제가 아니라, low-biomass 자료에서도 주 결론이 합리적인 전처리 범위에 의존하지 않음을 보이는 검증 자료로 제시하는 것이 적절하다.

Threshold 0.5는 control reads를 크게 줄이면서 50개 true sample과 91.64%의 true-sample reads를 유지한 균형점이다. Threshold 0.6 이상은 표본 손실이 커져 민감도 분석의 경계 조건으로만 제시해야 한다. 필터링이 강해질수록 공유 비율이 증가하는 현상은 sampling systems가 더 유사해졌다는 직접 증거가 아니라, 희귀 feature가 정의상 제거된 결과도 포함한다.

동일 50개 샘플의 matched-distance 분석은 이 주의점을 더 직접적으로 보여준다. 강한 필터는 Bray–Curtis와 unweighted UniFrac 기준의 apparent concordance를 높였지만 weighted UniFrac에서는 반대 방향의 변화를 보였다. 논문에서는 강한 filtering을 sampling-system agreement를 개선하는 방법으로 권고하기보다, 필터 선택이 agreement의 크기와 방향을 바꿀 수 있음을 제시해야 한다. 본 자료에서는 decontam threshold 0.5와 완화된 low-count filtering이 다섯 sampling system 비교의 원래 구조를 가장 안정적으로 보존했다.

## 재현 경로

- 통합 코드: `06_Integrated_decontam_filter_sampling_system_analysis.Rmd`
- 통합 결과: `Integrated_sensitivity_output/`
- 전체 검증표: `Integrated_sensitivity_output/Tables/Table16_final_sanity_checks.csv`
- 핵심 통합표: `Table14_integrated_decontam_threshold_results.csv`, `Table15_integrated_filter_results.csv`
- 동일 50개 샘플 검증: `Table17_fixed50_matched_distance_units.csv`–`Table20_fixed50_signal_ratios.csv`
- 추가 그림: `Figure7_fixed50_matched_distance_change`, `Figure8_fixed50_signal_preservation_ratio`
- Aitchison pseudocount 검증: `Table21_Aitchison_pseudocount_sensitivity.csv`
- 재사용 가능한 객체 및 결과: `Integrated_sensitivity_output/RDS/analysis_objects.rds`, `integrated_results.rds`, `fixed50_distance_results.rds`
- 실행 환경: `Integrated_sensitivity_output/sessionInfo.txt`

기존 논문 그림 코드와 결과인 `05_Final_main_supplement_figures.Rmd` 및 `Final_main_supplement_figure_output/`은 유지했다. 이전 singleton/doubleton/tenton 파생 Rmd와 06–09번 분산 분석 코드 및 각 결과 폴더는 본 통합 분석으로 대체했다.
