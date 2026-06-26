param(
    [string]$MatlabCommand = "matlab"
)

$ErrorActionPreference = "Stop"
$codeDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$scripts = @(
    "s031_build_spatial_confound_maps",
    "s032_control_spatial_confounds",
    "s033_spatial_autocorr_block_null",
    "s034_g1_binned_cortical_contribution",
    "s035_hemispheric_symmetry_analysis",
    "s036_summarize_local_story_extensions"
)

foreach ($script in $scripts) {
    Write-Host "Running $script"
    & $MatlabCommand -batch "cd('$codeDir'); $script"
    if ($LASTEXITCODE -ne 0) {
        throw "$script failed with exit code $LASTEXITCODE"
    }
}

Write-Host "Round-1 story extensions complete."
