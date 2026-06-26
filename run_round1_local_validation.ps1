$ErrorActionPreference = 'Stop'
$code='F:\Demo_engageX\ridge_gm_to_wm\code'
$logs='F:\Demo_engageX\ridge_gm_to_wm\logs'
$stages=@(
 's022_infer_local_r2_fc',
 's027_summarize_local_extension',
 's028_validate_local_r2_fc_across_sessions',
 's029_validate_gradient_difference_cross_session',
 's030_summarize_local_validations'
)
foreach($stage in $stages){
 $log=Join-Path $logs ($stage+'.log')
 $cmd="cd('$code'); addpath(genpath('$code')); diary('$log'); try, $stage; catch ME, disp(getReport(ME,'extended','hyperlinks','off')); diary off; exit(1); end; diary off; exit(0);"
 & matlab -batch $cmd
 if($LASTEXITCODE -ne 0){throw "$stage failed; see $log"}
}
