# Ligand library preparation: SMILES list -> 3D PDBQT via OpenBabel.
# Resumable (existing outputs are skipped); obabel exit codes are checked and
# failures are appended to failures.csv. Requires obabel on PATH; for conda
# builds set BABEL_DATADIR or --gen3d silently yields all-zero coordinates.
#
# Usage:  powershell -File prep_ligands.ps1 -LigandList ligands.smi -OutDir .
param(
    [string]$LigandList = "ligands.smi",   # lines: "SMILES<TAB>L0001"
    [string]$OutDir = ".",
    [string]$Obabel = "obabel"
)
$ErrorActionPreference = "Stop"
if (-not (Test-Path -LiteralPath $LigandList)) {
    Write-Error "missing input: $LigandList"
    exit 3
}
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$failCsv = Join-Path $OutDir "failures.csv"
$tmp = Join-Path $OutDir "one.smi"
$warnLog = Join-Path $OutDir "gen3d_warn.txt"

$i = 0
$nFail = 0
Get-Content -LiteralPath $LigandList | ForEach-Object {
    $i++
    if ([string]::IsNullOrWhiteSpace($_)) { return }
    $out = Join-Path $OutDir ("lig_{0:d4}.pdbqt" -f $i)
    if (Test-Path -LiteralPath $out) { return }          # resumable
    $_ | Set-Content -Encoding ascii -Path $tmp
    & $Obabel $tmp -O $out -h --gen3d 2>> $warnLog | Out-Null
    if ((-not (Test-Path -LiteralPath $out)) -or ($LASTEXITCODE -ne 0)) {
        "lig_{0:d4},obabel,rc={1}" -f $i, $LASTEXITCODE |
            Add-Content -Encoding utf8 -Path $failCsv
        $nFail++
    }
}
"prepared: {0} | failed: {1}" -f (Get-ChildItem -Path (Join-Path $OutDir "lig_*.pdbqt")).Count, $nFail
if ($nFail -gt 0) { exit 4 }
exit 0
