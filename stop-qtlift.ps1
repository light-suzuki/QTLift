$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$PidFile = Join-Path $Root 'data\runtime\qtlift.pid'

if (-not (Test-Path -LiteralPath $PidFile)) {
    Write-Host 'QTLift PID file does not exist; nothing was stopped.'
    exit 0
}

$processId = [int](Get-Content -LiteralPath $PidFile -Raw)
$process = Get-Process -Id $processId -ErrorAction SilentlyContinue
if ($process) {
    $commandLine = (Get-CimInstance Win32_Process -Filter "ProcessId=$processId").CommandLine
    if ($commandLine -notmatch 'QTLift[\\/]run\.py') {
        throw "PID $processId does not belong to QTLift; refusing to stop it."
    }
    Stop-Process -Id $processId
    Write-Host "Stopped QTLift PID $processId."
}
Remove-Item -LiteralPath $PidFile -Force

