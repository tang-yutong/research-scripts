# Rescue pass: re-run every ligand that still lacks a pose (timed out or
# failed in the parallel pass), this time serially with all cores per ligand
# and a shorter timeout. Failures are appended to failures.csv.
#
# Usage:  powershell -File dock_rescue.ps1 [-WorkDir .]
param(
    [string]$WorkDir = ".",
    [string]$Receptor = "receptor.pdbqt",
    [double]$CenterX = 30.0, [double]$CenterY = 30.0, [double]$CenterZ = 30.0,   # example
    [double]$SizeX = 24.0, [double]$SizeY = 24.0, [double]$SizeZ = 24.0,         # example
    [int]$Exhaustiveness = 8,
    [int]$NumModes = 9,
    [int]$TimeoutSec = 900,
    [int]$Cpu = 8,
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

Get-ChildItem -Path (Join-Path $WorkDir "lig_*.pdbqt") | ForEach-Object {
    $name = $_.BaseName
    $outPdbqt = Join-Path $WorkDir "out\$name.pdbqt"
    if (Test-Path -LiteralPath $outPdbqt) { return }
    $p = Start-Process -FilePath $Vina -ArgumentList @(
        '--receptor', $receptorPath, '--ligand', $_.FullName,
        '--center_x', "$CenterX", '--center_y', "$CenterY", '--center_z', "$CenterZ",
        '--size_x', "$SizeX", '--size_y', "$SizeY", '--size_z', "$SizeZ",
        '--cpu', "$Cpu", '--exhaustiveness', "$Exhaustiveness", '--num_modes', "$NumModes",
        '--out', $outPdbqt
    ) -RedirectStandardOutput (Join-Path $WorkDir "log\$name.log") `
      -RedirectStandardError (Join-Path $WorkDir "log\$name.err") `
      -PassThru -WindowStyle Hidden
    if (-not $p.WaitForExit($TimeoutSec * 1000)) {
        $p.Kill()
        "$name,vina-rescue,timeout" | Add-Content -Encoding utf8 -Path $failCsv
        if (Test-Path -LiteralPath $outPdbqt) { Remove-Item -LiteralPath $outPdbqt }
    } elseif ($p.ExitCode -ne 0) {
        "$name,vina-rescue,rc=$($p.ExitCode)" | Add-Content -Encoding utf8 -Path $failCsv
    }
}
exit 0
