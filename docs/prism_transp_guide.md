# PRISM TRANSP Viewer — Technical Guide

**PRISM v2.4.1** | 2026-04-14 | Jekil Lee (jklee@kfe.re.kr)

---

## Overview

`prism transp` (또는 `prism -t`)는 KSTAR 플라즈마 수송 해석 결과를 시각화하는 전용 뷰어입니다.
두 개의 독립적인 탭 그룹으로 구성됩니다:

| 그룹 | 탭 | 데이터 소스 | 설명 |
|------|-----|------------|------|
| **BiProfile** | Ti/vT Profile | MDS+ `biprofile` tree | BIPROFILE Bayesian inference 피팅 결과 |
| | ne/Te Profile | MDS+ `biprofile` tree | Thomson ne scale 적용, CES/Thomson raw overlay |
| | Ti/vT Time Trace | MDS+ `biprofile` tree | 선택된 ψ_N 위치에서의 시간 변화 |
| | ne/Te Time Trace | MDS+ `biprofile` tree | 선택된 ψ_N 위치에서의 시간 변화 |
| **TRANSP** | Profile | CDF 파일 (netCDF) | TRANSP OUTPUT CDF의 2D 프로파일 변수 |
| | Time Trace | CDF 파일 (netCDF) | TRANSP OUTPUT CDF의 1D 시계열 변수 |

---

## 실행 방법

```bash
# 다음 중 하나로 실행
prism -t
prism --transp
prism transp
```

---

## 사이드바 구조

```
 ▼ BiProfile
   ▼ Profiles
       Ti, vT
       ne, Te
   ▼ Time Traces
       Ti, vT
       ne, Te
 ▼ TRANSP
     Profile
     Time Trace
```

- 각 그룹/카테고리를 클릭하여 접기/펴기 가능 (▼/▶)
- 접힘 상태는 `~/.config/prism/settings.json`에 저장

---

## BiProfile 탭

### 데이터 소스

MDS+ 서버 (`mdsr.kstar.kfe.re.kr:8005`)의 `biprofile` tree에서 로드합니다.

| MDS+ Tree | 파라미터 |
|-----------|---------|
| `BITI` | Ion temperature (Ti) [keV] |
| `BIVT` | Toroidal rotation velocity (vT) [km/s] |
| `BITE` | Electron temperature (Te) [keV] |
| `BINE` | Electron density (ne) [10^19/m^3] |

각 파라미터에 대해 다음 노드를 읽습니다:
- `PSIN` — 정규화 폴로이달 자속 좌표
- `TIME` — 시간 배열 [s]
- `EFIT_USED` — 사용된 EFIT tree (e.g., EFIT01)
- `FIT_FUNC` — 피팅 함수명
- `MEAN`, `UNC` — 채널별 피팅 평균값/불확실성

### BiProfile Profile 탭 워크플로우

1. **Shot 입력** → Fetch 클릭
2. BIPROFILE + CES/Thomson raw + EFIT + DIAG_PARAMS를 병렬로 로드
3. **Available 리스트**에 time point 목록 표시 (e.g., `039551_005020 (Bi)`)
4. **Preview** 버튼으로 슬라이더 기반 데이터 탐색 (시간 이동, 재생, Fix Axes)
5. 원하는 time point를 **Selected**로 이동
6. **Plot** → ψ_N 기반 프로파일 표시
   - 피팅 결과: 실선 + 불확실성 밴드
   - Raw 데이터 overlay: CES (●), Thomson (■)
   - 사용/미사용 채널 구분 (미사용 = 회색)

### ne/Te Profile — TS ne Scale

- Thomson ne에 절대 보정 스케일 (scale factor)이 **기본 적용**됩니다
- "Unapply TS ne Scale" 토글로 스케일 해제 가능
- Browse Dialog에서도 동일 토글 제공

### BiProfile Time Trace 탭 워크플로우

1. Shot 입력 → Fetch
2. **Available 리스트**에 ψ_N 인덱스 목록 표시 (e.g., `039551_042 (Bi)`)
3. 원하는 ψ_N 위치를 Selected로 이동
4. Plot → 선택된 ψ_N 위치들의 시간 변화 표시 (평균 ± 불확실성)

### 플롯 옵션

| 옵션 | 기본값 | 설명 |
|------|--------|------|
| Color | Gradient(viridis) | 컬러맵 (Gradient 또는 Fixed) |
| Label font size | 12 | 축 라벨 글꼴 크기 |
| Legend font size | 8 | 범례 글꼴 크기 |
| Tick font size | 10 | 눈금 글꼴 크기 |
| Show Nodes | OFF | 피팅 노드 포인트 표시 |

---

## TRANSP CDF 탭

### 데이터 소스

TRANSP OUTPUT CDF 파일 (netCDF 포맷)을 로드합니다. 변수는 차원에 따라 자동 분류됩니다:

| 차원 | 분류 | 예시 |
|------|------|------|
| `(TIME,)` | Time Trace (1D) | Q0, PINJ, BETAT, PCUR, ... |
| `(TIME3, X)` | Profile (2D) | TE, TI, NE, NI, Q, CUR, CONDE, ... |
| `(TIME3, XB)` | Profile (2D) | 경계 정규화 좌표 기반 변수 |

**지원 포맷**: `.CDF`, `.cdf`, `.nc`
**netCDF 백엔드**: netCDF4 (기본) → scipy.io.netcdf (fallback)

### CDF 파일명 규칙

TRANSP CDF 파일명은 `{shot}{run_code}.CDF` 형태입니다:

| 예시 파일명 | Shot | Run Code |
|------------|------|----------|
| `39551W09.CDF` | 39551 | W09 |
| `41428A08.CDF` | 41428 | A08 |

### 기본 파일 경로

| 서버 | 기본 경로 |
|------|----------|
| nkstar | `/home/users/{USER}/` |
| ukstar (transp 그룹) | `/UKSTAR_HOME/data/transp/{USER}/` |
| ukstar (일반) | `/UKSTAR_HOME/{USER}/` |

### TRANSP Profile 탭 워크플로우

1. **Open CDF** → CDF 파일 선택 (non-modal 파일 다이얼로그)
2. Select Data 그룹 활성화, 상태 라벨에 `#39551 (W09) — 142 profiles, 85 time traces` 표시
3. **Filter**에 키워드 입력 → Variable 드롭다운 필터링
4. **Variable** 드롭다운에서 프로파일 변수 선택 (e.g., `TE - ELECTRON TEMPERATURE [KEV]`)
5. **Available 리스트**에 해당 변수의 time point 표시 (e.g., `39551_001500 (W09)`)
6. 원하는 time point를 **Selected**로 이동
7. **Plot** → 반경 방향 프로파일 표시

**다중 CDF 비교**:
- 다른 CDF를 Open하면 Available은 새 CDF의 time point로 갱신
- **Selected 리스트는 유지** — 이전 CDF에서 선택한 항목이 남아있음
- Plot 시 양쪽 CDF의 프로파일이 함께 표시됨
- 동일 CDF를 중복 로드하면 자동으로 무시

**범례 형식**: `#39551 001500ms (W09)`

### TRANSP Time Trace 탭 워크플로우

1. **Open CDF** → CDF 파일 선택
2. **Selected Runs**에 run ID 자동 추가 (e.g., `#39551 (W09)`)
3. **Filter** → Variable 드롭다운 필터링
4. **Variable** 선택 (e.g., `Q0 - SAFETY FACTOR ON AXIS`)
5. **Plot** → 선택된 run들의 해당 변수를 시간축으로 비교 표시

**다중 run 비교**:
- Open CDF를 반복하여 여러 run 추가 가능
- "Remove Selected Run" 버튼 또는 Delete 키로 run 제거 (CDF 캐시도 함께 해제)
- 동일 CDF는 중복 로드 불가

**범례 형식**: `#39551 (W09)`

### TRANSP 주요 변수 목록 (참고)

#### Profile 변수 (2D: time × radius)

| 변수명 | 설명 | 단위 |
|--------|------|------|
| TE | Electron temperature | keV |
| TI | Ion temperature | keV |
| NE | Electron density | m^-3 |
| NI | Ion density | m^-3 |
| Q | Safety factor | — |
| CUR | Total current density | A/cm^2 |
| CURB | Beam-driven current density | A/cm^2 |
| CURBS | Bootstrap current density | A/cm^2 |
| CONDE | Electron heat conductivity (χ_e) | — |
| CONDI | Ion heat conductivity (χ_i) | — |
| PLFLX | Poloidal flux | cm |

#### Time Trace 변수 (1D: time)

| 변수명 | 설명 | 단위 |
|--------|------|------|
| Q0 | Safety factor on axis | — |
| PINJ | Injected NBI power | W |
| BETAT | Toroidal beta | — |
| PCUR | Plasma current | A |
| TAUEE | Energy confinement time | s |
| NEUTT | Total neutron rate | /s |
| POHT | Ohmic heating power | W |
| VSURC | Surface voltage | V |

> 실제 변수 목록은 CDF 파일에 따라 다릅니다. Filter 기능으로 검색하세요.

---

## Settings 저장

모든 탭의 설정이 `~/.config/prism/settings.json`에 자동 저장됩니다:

| 탭 | Settings Key | 저장 항목 |
|----|-------------|----------|
| BiProfile Ti/vT Profile | `bi_tivt_profile` | shot, color_mode, font sizes, show_nodes, apply_scale |
| BiProfile ne/Te Profile | `bi_nete_profile` | shot, color_mode, font sizes, show_nodes, apply_scale |
| BiProfile Ti/vT Time Trace | `bi_tivt_timetrace` | shot, color_mode, font sizes, show_nodes |
| BiProfile ne/Te Time Trace | `bi_nete_timetrace` | shot, color_mode, font sizes, show_nodes |
| TRANSP Profile | `transp_profile` | color_mode, font sizes |
| TRANSP Time Trace | `transp_timetrace` | color_mode, font sizes |

---

## 키보드 단축키

| 키 | 동작 |
|----|------|
| Enter | Shot 입력 필드에서 Fetch 실행 (BiProfile) |
| Delete | Selected 리스트에서 선택 항목 제거 |

---

## 서버 환경

| 서버 | PRISM 경로 | Python |
|------|-----------|--------|
| nkstar | `/home/users/jklee/PRISM` | 3.8 (`/usr/bin/python38`) |
| ukstar | `/UKSTAR_HOME/jklee/PRISM` | 3.8 (`/usr/bin/python3.8`) |

**MDS+ 서버**: `mdsr.kstar.kfe.re.kr:8005` (내부망 전용)

---

## 파일 구조

```
PRISM/
├── data_loaders/
│   ├── biprofile_loader.py        # BiProfile MDS+ 로더
│   └── transp_cdf_loader.py       # TRANSP CDF (netCDF) 로더
├── ui/
│   ├── biprofile_profile_tab.py   # BiProfile Profile 탭
│   ├── biprofile_timetrace_tab.py # BiProfile Time Trace 탭
│   ├── transp_profile_tab.py      # TRANSP CDF Profile 탭
│   ├── transp_timetrace_tab.py    # TRANSP CDF Time Trace 탭
│   └── widgets/
│       └── biprofile_browse_dialog.py  # BiProfile Browse 다이얼로그
└── config/
    └── user_settings.py           # 사용자 설정 저장/로드
```

---

## 문의

- **개발자**: Jekil Lee (jklee@kfe.re.kr)
- **GitHub**: https://github.com/jekillee/PRISM
- **버전**: v2.4.1 (2026-04-14)
