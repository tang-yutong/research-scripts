# Parallel AutoDock Vina docking: N workers, each taking ligands whose index
# satisfies (index % N == W) -- start one process per W in 0..N-1. Per-ligand
# timeout kills hung jobs and quarantines the half-written pose; existing
# outputs are skipped (resumable); non-zero exits and timeouts land in
# failures.csv. Box centre/size below are EXAMPLES -- set them to your site.
#
# Usage:  powershell -File dock_parallel.ps1 -w 0 -n 8 [-WorkDir .]
param(
    [int]$w = 0,
    [int]$n = 8,
    [string]$WorkDir = ".",                # holds lig_*.pdbqt; writes out/, log/
    [string]$Receptor = "receptor.pdbqt",
    [double]$CenterX = 30.0, [double]$CenterY = 30.0, [double]$CenterZ = 30.0,   # example
    [double]$SizeX = 24.0, [double]$SizeY = 24.0, [double]$SizeZ = 24.0,         # example
    [int]$Exhaustiveness = 8,
    [int]$NumModes = 9,
    [int]$TimeoutSec = 1200,
    [string]$Vina = "vina"
)
$ErrorActionPreference = "Stop"
$receptorPath = Join-Path $WorkDir $Receptor
if (-not (Test-Path -LiteralPath $receptorPath)) {
    Write-Error "missing receptor: $receptorPath"
    exit 3
}
New-Item -ItemType Directory -Force -Path (Join-Path $WorkDir "out"),
    (Join-Path $WorkDir "log") | Out-Null
$failCsv = Join-Path $WorkDir "failures.csv"

$i = -1
Get-ChildItem -Path (Join-Path $WorkDir "lig_*.pdbqt") | ForEach-Object {
    $i++
    if (($i % $n) -ne $w) { return }
    $name = $_.BaseName
    $outPdbqt = Join-Path $WorkDir "out\$name.pdbqt"
    if (Test-Path -LiteralPath $outPdbqt) { return }     # resumable
    $p = Start-Process -FilePath $Vina -ArgumentList @(
        '--receptor', $receptorPath, '--ligand', $_.FullName,
        '--center_x', "$CenterX", '--center_y', "$CenterY", '--center_z', "$CenterZ",
        '--size_x', "$SizeX", '--size_y', "$SizeY", '--size_z', "$SizeZ",
        '--cpu', '1', '--exhaustiveness', "$Exhaustiveness", '--num_modes', "$NumModes",
        '--out', $outPdbqt
    ) -RedirectStandardOutput (Join-Path $WorkDir "log\$name.log") `
      -RedirectStandardError (Join-Path $WorkDir "log\$name.err") `
      -PassThru -WindowStyle Hidden
    if (-not $p.WaitForExit($TimeoutSec * 1000)) {
        $p.Kill()
        "$name,vina,timeout" | Add-Content -Encoding utf8 -Path $failCsv
        if (Test-Path -LiteralPath $outPdbqt) { Remove-Item -LiteralPath $outPdbqt }
    } elseif ($p.ExitCode -ne 0) {
        "$name,vina,rc=$($p.ExitCode)" | Add-Content -Encoding utf8 -Path $failCsv
    }
}
exit 0
