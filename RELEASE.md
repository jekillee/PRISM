# PRISM Release Guide

버전 릴리즈 및 배포 절차 안내

## 버전 업데이트 체크리스트

새 버전 릴리즈 시 다음 파일들을 업데이트:

1. `config/app_config.py` - VERSION 변수
2. `CHANGELOG.md` - 변경 내역 추가

## 배포 절차

### 1. 코드 수정 및 테스트

인터넷망 PC에서 코드 수정 후 테스트 완료

### 2. 연구원 내부망 서버 배포

내부망은 GitHub 연결이 안 되므로 파일을 직접 복사:

```bash
# 방법 A: USB로 전체 폴더 복사
# 인터넷망 -> USB -> 내부망

# 방법 B: 변경된 파일만 복사 (권장)
# 변경 파일 확인 후 해당 파일만 복사
```

내부망 서버 경로: `/home/users/jklee/PRISM`

### 3. 내부망 Git 커밋

내부망 터미널 (nkstar)에서:

```bash
cd /home/users/jklee/PRISM

# 변경사항 확인
git status

# 스테이징
git add -A

# 커밋
git commit -m "v1.1.2: 변경 내용 요약"

# 태그 추가
git tag -a v1.1.2 -m "v1.1.2"
```

**참고**: 내부망에서는 `git push`가 안 됨 (GitHub 연결 불가)

### 4. GitHub 배포

인터넷망 PC (PowerShell)에서:

```powershell
cd "D:\Research\Code\Code Scripts\PRISM\v1.1.2"

# 변경사항 확인
git status

# 스테이징 및 커밋
git add -A
git commit -m "v1.1.2: 변경 내용 요약"

# 태그 추가
git tag -a v1.1.2 -m "v1.1.2"

# GitHub에 push
git push origin master
git push origin v1.1.2
```

## 서버 정보

| 환경 | 경로/주소 |
|------|-----------|
| 내부망 서버 | `/home/users/jklee/PRISM` |
| GitHub | `https://github.com/jekillee/PRISM.git` |
| Python | `/usr/bin/python38` (내부망) |

## 커밋 메시지 형식

```
v{버전}: 변경 내용 요약

예시:
v1.1.2: Fix MSE Time Trace gamma_min error and q/j fill_between style
v1.1.1: UI improvements for N-Mode Spectrum, TV, TivT tabs
v1.1.0: Add user settings persistence and update notification
```
