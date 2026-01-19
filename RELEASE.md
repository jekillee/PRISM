# PRISM Release Guide

Version release and deployment procedures.

## Version Update Checklist

Update the following files for each new release:

1. `config/app_config.py` - VERSION variable
2. `CHANGELOG.md` - Add change history

## Deployment Procedure

### 1. Code Modification and Testing

Modify and test code on the internet-connected PC.

### 2. Deploy to KFE Internal Server

The internal network cannot connect to GitHub, so files must be copied manually:

```bash
# Method A: Copy entire folder via USB
# Internet PC -> USB -> Internal network

# Method B: Copy only modified files (recommended)
# Check modified files and copy only those
```

Internal server path: `/home/users/jklee/PRISM`

### 3. Git Commit on Internal Network

On the internal terminal (nkstar):

```bash
cd /home/users/jklee/PRISM

# Check changes
git status

# Stage changes
git add -A

# Commit
git commit -m "v1.1.2: Brief description of changes"

# Add tag
git tag -a v1.1.2 -m "v1.1.2"
```

**Note**: `git push` does not work on the internal network (no GitHub connection).

### 4. Deploy to GitHub

On the internet-connected PC (PowerShell):

```powershell
cd "D:\Research\Code\Code Scripts\PRISM\v1.1.2"

# Check changes
git status

# Stage and commit
git add -A
git commit -m "v1.1.2: Brief description of changes"

# Add tag
git tag -a v1.1.2 -m "v1.1.2"

# Push to GitHub
git push origin master
git push origin v1.1.2
```

## Server Information

| Environment | Path/URL |
|-------------|----------|
| Internal Server | `/home/users/jklee/PRISM` |
| GitHub | `https://github.com/jekillee/PRISM.git` |
| Python | `/usr/bin/python38` (internal) |

## Commit Message Format

```
v{version}: Brief description of changes

Examples:
v1.1.2: Fix MSE Time Trace gamma_min error and q/j fill_between style
v1.1.1: UI improvements for N-Mode Spectrum, TV, TivT tabs
v1.1.0: Add user settings persistence and update notification
```
