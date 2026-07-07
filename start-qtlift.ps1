param(
    [switch]$Rebuild,
    [switch]$NoBrowser
)

$ErrorActionPreference = 'Stop'
$Root = Split-Path -Parent $MyInvocation.MyCommand.Path
$Python = Join-Path $Root '.venv\Scripts\python.exe'
$Frontend = Join-Path $Root 'frontend'
$DistIndex = Join-Path $Frontend 'dist\index.html'
$Runtime = Join-Path $Root 'data\runtime'
$PidFile = Join-Path $Runtime 'qtlift.pid'
$Url = 'http://127.0.0.1:8765'

New-Item -ItemType Directory -Force -Path $Runtime | Out-Null

function Test-QTLiftHealth {
    try {
        $health = Invoke-RestMethod -Uri "$Url/api/health" -TimeoutSec 2
        return $health.status -eq 'ok'
    } catch {
        return $false
    }
}

function Open-QTLiftBrowser {
    if ($NoBrowser) { return }
    $chromeCandidates = @(
        "$env:ProgramFiles\Google\Chrome\Application\chrome.exe",
        "${env:ProgramFiles(x86)}\Google\Chrome\Application\chrome.exe",
        "$env:LOCALAPPDATA\Google\Chrome\Application\chrome.exe"
    )
    $chrome = $chromeCandidates | Where-Object { $_ -and (Test-Path -LiteralPath $_) } | Select-Object -First 1
    if ($chrome) {
        Start-Process -FilePath $chrome -ArgumentList "--app=$Url" | Out-Null
    } else {
        Start-Process $Url | Out-Null
    }
}

if (Test-QTLiftHealth) {
    Write-Host "QTLift is already running at $Url"
    Open-QTLiftBrowser
    exit 0
}

if (-not (Test-Path -LiteralPath $Python)) {
    py -3.13 -m venv (Join-Path $Root '.venv')
    if ($LASTEXITCODE -ne 0) { throw 'Failed to create the QTLift virtual environment.' }
    & $Python -m pip install -r (Join-Path $Root 'requirements.txt')
    if ($LASTEXITCODE -ne 0) { throw 'Failed to install QTLift Python dependencies.' }
}

$sourceFiles = Get-ChildItem -LiteralPath (Join-Path $Frontend 'src') -Recurse -File
$needsBuild = $Rebuild -or -not (Test-Path -LiteralPath $DistIndex)
if (-not $needsBuild) {
    $distTime = (Get-Item -LiteralPath $DistIndex).LastWriteTimeUtc
    $needsBuild = [bool]($sourceFiles | Where-Object { $_.LastWriteTimeUtc -gt $distTime } | Select-Object -First 1)
}
if ($needsBuild) {
    Push-Location $Frontend
    try {
        if (-not (Test-Path -LiteralPath (Join-Path $Frontend 'node_modules'))) {
            npm install
            if ($LASTEXITCODE -ne 0) { throw 'npm install failed.' }
        }
        npm run build
        if ($LASTEXITCODE -ne 0) { throw 'Frontend build failed.' }
    } finally {
        Pop-Location
    }
}

$stdout = Join-Path $Runtime 'server.stdout.log'
$stderr = Join-Path $Runtime 'server.stderr.log'
$process = Start-Process -FilePath $Python -ArgumentList (Join-Path $Root 'run.py') `
    -WorkingDirectory $Root -WindowStyle Hidden -PassThru `
    -RedirectStandardOutput $stdout -RedirectStandardError $stderr
for ($attempt = 0; $attempt -lt 30; $attempt++) {
    Start-Sleep -Milliseconds 500
    if (Test-QTLiftHealth) {
        $listenerLine = netstat -ano | Select-String '127\.0\.0\.1:8765\s+.*LISTENING\s+\d+' | Select-Object -First 1
        if (-not $listenerLine) { throw 'QTLift is healthy but its listener PID could not be resolved.' }
        $listenerPid = [int](($listenerLine.Line -split '\s+')[-1])
        Set-Content -LiteralPath $PidFile -Value $listenerPid -Encoding ascii
        Write-Host "QTLift started at $Url (PID $listenerPid)"
        Open-QTLiftBrowser
        exit 0
    }
    if ($process.HasExited) { break }
}

$errorTail = if (Test-Path -LiteralPath $stderr) { Get-Content -LiteralPath $stderr -Tail 30 | Out-String } else { '' }
throw "QTLift failed to start.`n$errorTail"
