# CSP + Stable Lattice + XRD Pipeline

입력 구조(CIF/POSCAR 등)를 넣으면 아래를 자동으로 수행하는 재현 가능한 Python CLI 프로젝트입니다.

1. **CSP-like 후보 생성/스코어링**: 여러 브라베 격자 템플릿(정방정계, 육방정계 등)으로 구조 후보를 만들고 pseudo-energy 점수를 계산
2. **가장 안정한 결정형 격자 선택**: 최저 score 후보를 best structure로 선택
3. **XRD 시뮬레이션**: 선택된 구조의 powder XRD 패턴(Cu Kα)을 CSV/PNG로 저장

> 주의: 이 코드는 "빠른 스크리닝/프로토타입" 목적입니다. DFT 기반 정밀 CSP를 대체하지 않습니다.

## 빠른 시작

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

실행:

```bash
cspxrd --input examples/nacl.cif --output-dir outputs/nacl
```

생성 결과:

- `outputs/nacl/report.json`
- `outputs/nacl/predicted_best_structure.cif`
- `outputs/nacl/xrd_pattern.csv`
- `outputs/nacl/xrd_pattern.png`

## GitHub에서 재생 가능하게 쓰는 방법

1. 이 폴더를 GitHub repo에 push
2. 아래 명령만으로 같은 결과 재현

```bash
git clone <your-repo-url>
cd <repo>
python -m venv .venv
source .venv/bin/activate
pip install -e .
cspxrd --input examples/nacl.cif --output-dir outputs/nacl
```

## 테스트

```bash
pip install -e .[dev]
pytest
```
