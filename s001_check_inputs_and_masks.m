%% s001_check_inputs_and_masks
clear; clc;
cfg = load_project_config();
out_dir = cfg.result_dirs{2};

atlas_info = niftiinfo(cfg.schaefer400_file);
wm_info = niftiinfo(cfg.wm_mask_source);
gm_info = niftiinfo(cfg.gm_mask_source);
atlas = single(niftiread(atlas_info));
wm = niftiread(wm_info) > 0;
gm = niftiread(gm_info) > 0;

[atlas_wm_ok, atlas_wm_msg] = assert_same_space(atlas_info, wm_info);
[gm_wm_ok, gm_wm_msg] = assert_same_space(gm_info, wm_info);
if ~atlas_wm_ok || ~gm_wm_ok
    error('Mask/atlas space mismatch: atlas-WM %s; GM-WM %s', ...
        atlas_wm_msg, gm_wm_msg);
end

if cfg.erode_wm_mask
    if ~license('test','image_toolbox')
        error('WM erosion requires Image Processing Toolbox.');
    end
    wm = imerode(wm, strel('sphere', cfg.wm_erosion_radius_vox));
end

atlas_ids = unique(atlas(atlas > 0));
missing_roi = setdiff((1:400)', double(atlas_ids(:)));
overlap_atlas_wm = nnz((atlas > 0) & wm);
overlap_gm_wm = nnz(gm & wm);

% Write the exact mask used by all later stages.
out_base = erase(cfg.wm_mask_final, '.nii.gz');
niftiwrite(single(wm), out_base, wm_info, 'Compressed', true);

rows = {};
if isfield(cfg,'all_subject_list')
    qc_subject_list = cfg.all_subject_list;
else
    qc_subject_list = cfg.subject_list;
end
for s = 1:numel(qc_subject_list)
    subject_id = qc_subject_list{s};
    for r = 1:numel(cfg.run_dirs)
        f = get_run_file(cfg, subject_id, r);
        exists_flag = exist(f,'file') == 2;
        n_time = NaN; tr = NaN; same_space = NaN;
        header_checked = s <= cfg.qc_header_subjects_per_run;
        status = "missing";
        if exists_flag
            if header_checked
                info = niftiinfo(f);
                n_time = info.ImageSize(4);
                tr = info.PixelDimensions(4);
                [same_space,~] = assert_same_space(info, wm_info);
                if same_space && n_time == cfg.expected_timepoints && ...
                        abs(tr-cfg.TR) < 1e-3
                    status = "valid_checked";
                else
                    status = "mismatch";
                end
            else
                % Every run is fully checked again when s002 reads it.
                status = "exists_pending_s002_header_check";
            end
        end
        rows(end+1,:) = {subject_id,cfg.run_ids{r},f,exists_flag, ...
            header_checked,n_time,tr,same_space,status}; %#ok<SAGROW>
    end
end
Q = cell2table(rows, 'VariableNames', {'subject_id','run_id','file', ...
    'exists','header_checked','n_timepoints','TR','same_space','status'});
writetable(Q, fullfile(out_dir,'valid_subject_run_table.csv'));

metric = ["atlas_wm_space";"gm_wm_space";"wm_voxels";"gm_voxels"; ...
    "atlas_voxels";"atlas_wm_overlap";"gm_wm_overlap";"missing_rois"; ...
    "existing_runs"];
value = [string(atlas_wm_msg);string(gm_wm_msg);string(nnz(wm)); ...
    string(nnz(gm));string(nnz(atlas>0));string(overlap_atlas_wm); ...
    string(overlap_gm_wm);strjoin(string(missing_roi),';'); ...
    string(nnz(Q.exists))];
summary = table(metric,value);
writetable(summary, fullfile(out_dir,'QC_mask_alignment_report.csv'));

fprintf('WM voxels: %d\n', nnz(wm));
fprintf('Atlas-WM overlap: %d; GM-WM overlap: %d\n', ...
    overlap_atlas_wm, overlap_gm_wm);
fprintf('Existing runs: %d/%d; headers sampled: %d\n', ...
    nnz(Q.exists),height(Q),nnz(Q.header_checked));
if ~isempty(missing_roi)
    error('Atlas is missing ROI labels: %s', strjoin(string(missing_roi),','));
end
if any(~Q.exists) || any(Q.status == "mismatch")
    warning('Some runs failed QC. Inspect valid_subject_run_table.csv.');
end
