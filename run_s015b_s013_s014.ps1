$ErrorActionPreference = 'Stop'

$project = 'F:\Demo_engageX\ridge_gm_to_wm'
$code = Join-Path $project 'code'
$logDir = Join-Path $project 'logs'
$masterLog = Join-Path $logDir 's015b_s013_s014_master.log'

if (-not (Test-Path -LiteralPath $logDir)) {
    New-Item -ItemType Directory -Path $logDir | Out-Null
}

$stages = @(
    's015b_compare_beta_fc_gradient_stability',
    's013_run_lagged_ridge_model',
    's014_null_model_ridge_prediction'
)

"Pipeline started: $(Get-Date -Format o)" | Set-Content -LiteralPath $masterLog

foreach ($stage in $stages) {
    $stageLog = Join-Path $logDir ($stage + '.rerun.log')
    "START $stage $(Get-Date -Format o)" | Add-Content -LiteralPath $masterLog
    $matlabCommand = "cd('$code'); addpath(genpath('$code')); diary('$stageLog'); try, $stage; catch ME, disp(getReport(ME,'extended','hyperlinks','off')); diary off; exit(1); end; diary off; exit(0);"
    & matlab -batch $matlabCommand
    $exitCode = $LASTEXITCODE
    "END $stage exit=$exitCode $(Get-Date -Format o)" | Add-Content -LiteralPath $masterLog
    if ($exitCode -ne 0) {
        throw "Stage $stage failed with exit code $exitCode. See $stageLog"
    }
}

"Pipeline completed: $(Get-Date -Format o)" | Add-Content -LiteralPath $masterLog
