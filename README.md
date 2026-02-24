# CSP + Stable Lattice + XRD Pipeline

이 프로젝트는 아래 2가지 실행 모드를 제공합니다.

1. **CLI 모드**: 구조 파일(CIF/POSCAR)을 입력해 CSP-like 스코어링 + XRD 생성
2. **웹 모드**: **SMILES 문자열 입력**으로 구조를 생성하고 XRD 그래프를 웹에서 바로 확인

> 주의: 이 코드는 빠른 프로토타입/스크리닝 목적이며, DFT 기반 정밀 CSP를 대체하지 않습니다.

## 설치

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

## 1) CLI 사용법

```bash
cspxrd --input examples/nacl.cif --output-dir outputs/nacl
```

생성 결과:
- `report.json`
- `predicted_best_structure.cif`
- `xrd_pattern.csv`
- `xrd_pattern.png`

## 2) 웹 사용법 (SMILES → XRD)

```bash
cspxrd-web
```

브라우저에서 `http://localhost:8000` 접속 후 SMILES 입력 (예: `CCO`, `c1ccccc1`).

## GitHub에서 재생

```bash
git clone <your-repo-url>
cd <repo>
python -m venv .venv
source .venv/bin/activate
pip install -e .
cspxrd-web
```

## 테스트

```bash
python -m compileall src tests
# 의존성 설치된 환경에서
pytest
```
