$ErrorActionPreference = 'Stop'

$project = 'F:\Demo_engageX\ridge_gm_to_wm'
$code = Join-Path $project 'code'
$logDir = Join-Path $project 'logs'
if (-not (Test-Path -LiteralPath $logDir)) {
    New-Item -ItemType Directory -Path $logDir | Out-Null
}

$stages = @(
    's004b_refit_shared_lambda_betas',
    's005_write_r2_maps',
    's006_extract_beta_profiles',
    's007_write_roi_beta_maps',
    's008_compute_network_contribution',
    's010_compute_network_specificity',
    's011_compute_fc_baseline',
    's011b_compare_fc_vs_beta',
    's012_run_network_level_ridge_model',
    's009b_beta_kmeans_parcellation',
    's009a_brainspace_gradient',
    's015a_compare_beta_fc_stability',
    's017_behavior_prediction',
    's015b_compare_beta_fc_gradient_stability'
)

$masterLog = Join-Path $logDir 'remaining_pipeline_master.log'
"Pipeline started: $(Get-Date -Format o)" | Add-Content -LiteralPath $masterLog

foreach ($stage in $stages) {
    $stageLog = Join-Path $logDir ($stage + '.log')
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
