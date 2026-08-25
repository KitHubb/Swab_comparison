# DT-Swab 원고 참고문헌·벤치마크 병목 점검

- 점검일: 2026-08-25
- 대상 원고: `01_CoST_공저자검토용_통합원고_20260726_ksy2.docx`
- 원칙: 기존 sampling-system 비교라는 주제를 유지하고, low-biomass contamination 분석은 결론의 신뢰도를 보강하는 근거로 사용한다.

## 1. 현재 가장 큰 병목

분석 자체는 완료되어 있다. 현재 병목은 분석 실행이 아니라 다음 세 항목의 원고 동기화다.

1. 본문 인용번호 일부가 현재 References 목록과 맞지 않는다.
2. `decontam`과 low-biomass contamination의 핵심 방법론 문헌이 References에 없다.
3. 통계 및 전처리 문장이 최종 분석 설정과 일치하지 않는다.

이 상태에서는 sampling-system 비교 결과는 제시할 수 있지만, contamination-aware preprocessing을 연구의 강점으로 안전하게 주장하기 어렵다.

## 2. 기존 참고문헌과 현 연구의 비교

| 문헌 | 목적·데이터 | 주 분석 | 주요 결과 | 현 연구와의 관계 |
|---|---|---|---|---|
| Grice et al., 2008 | 건강인 antecubital fossa에서 swab, scrape, biopsy 비교 | 16S OTU 구성과 군집 비교 | 채취 깊이가 달라도 우점 군집은 대체로 유사했고 개인차가 큼 | 현 연구는 침습도가 다른 방법이 아니라 5개 완전한 swab–medium system을 비교함 |
| Bjerre et al., 2019 | 건강인 9명, eSwab/scrape 및 12개 추출 키트 | 16S, library 성공률, diversity, overlap | eSwab과 scrape sequence overlap 99.3%; 추출 키트의 영향이 큼 | 현 연구는 downstream extraction을 고정하여 sampling system 효과에 집중함 |
| Ogai et al., 2018 | 건강인 7명, swab 대 tape stripping | 16S와 배양 | 상대조성은 유사하지만 tape에서 배양 가능한 균이 더 많이 회수됨 | 현 연구는 상대조성 비교이며 viability·절대 회수율은 평가하지 않음 |
| Balacco et al., 2025 | 16명, cotton/eSwab, 습윤액, 시간, 보관조건 비교 | 16S, DNA yield, mixed model | eSwab total DNA가 높았지만 community profile 차이는 작고 개인차가 큼 | 가장 가까운 benchmark. 현 연구는 5개 system과 두 피부 부위를 평가하지만 n=5임 |
| Manus et al., 2022 | 영아 9명, 3개 부위, wet/dry와 ethanol/freezer 비교 | 16S, PERMANOVA, mixed-effects model | 처리조건 약 10%, 부위 약 20% 설명 | 현 연구는 즉시 냉동을 고정하고 system 효과를 평가함 |
| Panpradist et al., 2014 | 7개 swab의 세균 방출 효율 실험 | 단일 세균 qPCR | swab 재질·시료량·회수 조건에 따라 효율이 달라짐 | swab 물성의 간접 근거이며 피부 community 결과의 직접 근거는 아님 |
| Wise et al., 2021 | 4개 swab에서 단일 세균 DNA 회수 | qPCR | flocked swab의 회수량이 높음 | 통제된 forensic 실험으로, 피부 상대조성에 직접 외삽하면 안 됨 |
| Bruijns et al., 2018 | 5개 swab에서 순수 DNA 회수 | DNA extraction/recovery efficiency | 모든 swab에서 회수율이 50% 미만, nylon-flocked가 상대적으로 우수 | swab 물성 설명에만 제한적으로 사용 |
| Klymiuk et al., 2016 | 건강인 8명, 3개 부위, −80°C 장기 보관 | 16S diversity와 taxonomic ratio | richness/Shannon은 안정적이나 일부 분류군 비율은 변함 | 현 연구의 즉시 냉동 조건을 벗어난 운반·장기 안정성은 주장할 수 없음을 보여줌 |
| Kim et al., 2026 | 건강한 한국인 10명, 8개 피부 부위 | V1–V3 16S, multivariate analysis | 해부학적 부위 효과가 생리학적 skin type 효과보다 큼 | 현 연구의 participant/site 효과가 system 효과보다 크다는 결과를 뒷받침함 |

## 3. 현 연구의 방어 가능한 novelty

현 연구의 차별점은 새로운 swab 소재를 개발했다는 데 있지 않다. 다음 조합이 비교 연구로서의 기여다.

1. 실제 임상연구에서 사용할 수 있는 5개 integrated swab–medium system을 비교했다.
2. 동일 참가자에게 모든 system을 적용하고 sebaceous site와 moist site를 함께 평가했다.
3. 모든 시료에 동일한 저장, 추출, library preparation, sequencing 및 bioinformatic workflow를 적용했다.
4. low-biomass negative control과 `decontam` threshold 평가를 포함해 전처리 선택이 sampling-system 결론을 바꾸는지 확인했다.

권장 positioning:

> Under a standardized immediate-freezing and downstream processing workflow, the predominant bacterial profiles were comparatively stable across five integrated swab–medium systems, whereas participant and anatomical-site effects were larger. Contamination-aware and filtering-sensitivity analyses supported the stability of this interpretation in low-biomass skin samples.

피해야 할 표현:

- `the five systems are equivalent`
- `the systems are fully interchangeable`
- `a universal standardized sampling method was established`
- `the optimal swab was identified`

## 4. 인용번호 불일치 점검

| 원고 위치 | 현재 인용 | 문제 | 필요한 조치 |
|---|---|---|---|
| Introduction 첫 단락 | `[1,2]` | 현재 1, 2번은 QIIME 2와 DADA2이며 anatomical-site·inter-individual variation 근거가 아님 | Grice, Kim 및 적절한 skin microbiome review로 교체 |
| Introduction 첫 단락 | `[3,6,7]` | SILVA와 taxonomy classifier 문헌은 pre-analytical variation 근거가 아님 | Bjerre, Balacco, Manus 및 low-biomass contamination 문헌으로 교체 |
| Introduction 둘째 단락 | `[1,3,4]` | collection-method 비교 근거와 번호가 맞지 않음 | Grice, Bjerre, Ogai로 교체 |
| Introduction 셋째 단락 | `[5,8]` | 5번 q2-feature-classifier는 swabbing condition 근거가 아님 | Bjerre, Balacco, Manus를 사용 |
| Participants/Methods | `[2,6]` | DADA2와 taxonomy 논문은 pre-sampling restriction 근거가 아님 | 실제 적용한 임상 protocol 또는 Kim et al. 문헌으로 교체 |
| Discussion swab comfort | `[11]` | 11번 Manus 연구는 pediatric nasopharyngeal comfort 연구가 아님 | 해당 원문을 새로 추가하거나 문장 삭제 |
| Discussion eSwab viability | `[15]` | 15번 Klymiuk 연구는 −80°C 보관에 관한 피부 16S 연구이며 eSwab culture viability 근거가 아님 | Amies/eSwab viability 원문을 추가하거나 주장을 제한 |
| Negative controls | 없음 | `decontam` 방법을 사용하지만 원 논문 인용이 없음 | Davis et al., 2018 추가 |
| Low-biomass 배경·한계 | 없음 | contamination-aware analysis의 필요성을 뒷받침하는 문헌 없음 | Eisenhofer et al., 2019 등 추가 |

## 5. 분석 설정과 원고의 불일치

1. Sequence preprocessing에 `without fixed-length truncation (trunc-len-f = 270, trunc-len-r = 240)`라고 되어 있다. 270/240은 고정 truncation이므로 `with fixed truncation lengths of 270 and 240 bp`로 고쳐야 한다.
2. Statistical Analysis에는 Bonferroni adjustment로 적혀 있다. 최종 sensitivity 분석에서 BH FDR을 사용했다면 본 분석·표·legend와 함께 하나의 방식으로 통일해야 한다.
3. Taxonomic Assignment의 `A naïve Bayes classifier was trained using pre-trained the SILVA...` 문장은 trained classifier인지 pre-trained classifier인지 모순된다.
4. Results에서 `Figures S4–S7`을 인용하지만 현재 legend 목록은 Figure S4–S6까지만 확인된다.
5. `collection medium`의 대소문자 표기가 Methods 안에서 일관되지 않다.
6. `shared dominant taxa`라는 표현은 genus-level 분석이 abundance filtering 없이 수행되었다면 정확하지 않다. `shared genus-level taxa`가 더 안전하다.

## 6. References에 추가할 최소 문헌

1. Davis NM, Proctor DM, Holmes SP, Relman DA, Callahan BJ. Simple statistical identification and removal of contaminant sequences in marker-gene and metagenomics data. Microbiome. 2018;6:226. doi:10.1186/s40168-018-0605-2.
2. Eisenhofer R, Minich JJ, Marotz C, Cooper A, Knight R, Weyrich LS. Contamination in low microbial biomass microbiome studies: issues and recommendations. Trends Microbiol. 2019;27:105–117. doi:10.1016/j.tim.2018.11.003.

## 7. 다음 실행 순서

1. References에 Davis와 Eisenhofer 문헌 추가
2. 본문 citation numbering 전체 재매핑
3. DADA2 truncation, BH FDR 및 taxonomy classifier 문장 동기화
4. Figure/Table 번호와 legend 대조
5. 수정된 원고와 sensitivity report의 수치 교차검증
6. 최종 원고에서 sampling-system comparison을 주 결과로 유지하고 decontam/filter sensitivity는 신뢰도 보강 결과로 배치

