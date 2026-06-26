$ErrorActionPreference = 'Stop'
$code='F:\Demo_engageX\ridge_gm_to_wm\code'
$logs='F:\Demo_engageX\ridge_gm_to_wm\logs'
$stages=@(
 's100_config_round2',
 's101_qc_round2_masks',
 's102_extract_round2_timeseries',
 's103_build_round2_manifests',
 's104_run_round2_zerolag',
 's105_run_round2_lag10',
 's106_run_round2_lag_smooth_basis',
 's107_summarize_round2_lag_decision'
)
foreach($stage in $stages){
 $log=Join-Path $logs ($stage+'.log')
 $cmd="cd('$code'); addpath(genpath('$code')); diary('$log'); try, $stage; catch ME, disp(getReport(ME,'extended','hyperlinks','off')); diary off; exit(1); end; diary off; exit(0);"
 & matlab -batch $cmd
 if($LASTEXITCODE -ne 0){throw "$stage failed; see $log"}
}
