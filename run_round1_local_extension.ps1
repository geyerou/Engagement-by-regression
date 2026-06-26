$ErrorActionPreference = 'Stop'
$code='F:\Demo_engageX\ridge_gm_to_wm\code'
$logs='F:\Demo_engageX\ridge_gm_to_wm\logs'
$stages=@(
 's021_compute_local_r2_fc',
 's022_infer_local_r2_fc',
 's023_build_beta_fc_gradient_templates',
 's024_infer_subject_beta_fc_gradient_difference',
 's025_localize_local_metrics_jhu',
 's026_predict_behavior_from_local_metrics',
 's027_summarize_local_extension'
)
foreach($stage in $stages){
 $log=Join-Path $logs ($stage+'.log')
 $cmd="cd('$code'); addpath(genpath('$code')); diary('$log'); try, $stage; catch ME, disp(getReport(ME,'extended','hyperlinks','off')); diary off; exit(1); end; diary off; exit(0);"
 & matlab -batch $cmd
 if($LASTEXITCODE -ne 0){throw "$stage failed; see $log"}
}
