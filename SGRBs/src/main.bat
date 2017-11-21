del run.exe *.mod *.obj *.pdb *.ilk *.smod

REM ifort /debug /standard-semantics /traceback /gen-interfaces /check /fpe:0 main.f90 -o run.exe
ifort /O2 /standard-semantics /assume:realloc_lhs main.f90 -o run.exe
REM ifort /O3 /standard-semantics /assume:realloc_lhs main.f90 -o run.exe

@echo off
TIMEOUT 1 >nul 2>nul
editbin /STACK:99999999 run.exe
TIMEOUT 1 >nul 2>nul
@echo on
@echo off