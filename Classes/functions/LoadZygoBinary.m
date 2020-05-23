function data = LoadZygoBinary(filename)
% data = LoadZygoBinary(filename)
%
% loads all data from a Zygo binary (.dat) file to a structure; no rescaling is performed
% More technical data can be found in the Metropro Reference Guide
%
% Author: Massimo Galimberti 2010 (a pure awesome job!)
%
%%% open file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[fid, msg] = fopen(filename, 'r'); %#ok<NASGU>
if fid<0
    error('msg');
end

%
%%% HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

data.magic_number       = fread(fid, 1,  'int32',   0, 'b');
data.header_format      = fread(fid, 1,  'int16',   0, 'b');
data.header_size        = fread(fid, 1,  'int32',   0, 'b');
% /!\ This is for Durango, who gives a wrong header_size
% header size can be 4096 (newer versions) or 834 (older versions)
if data.header_size ~= 4096
    data.header_size = 834;
end
data.swinfo.type        = fread(fid, 1,  'int16',   0, 'b');
data.swinfo.date        = fread(fid, 30, '*char'         )';
data.swinfo.vers.maj    = fread(fid, 1,  'int16',   0, 'b');
data.swinfo.vers.min    = fread(fid, 1,  'int16',   0, 'b');
data.swinfo.vers.bug    = fread(fid, 1,  'int16',   0, 'b');
data.ac_org_x           = fread(fid, 1,  'int16',   0, 'b');
data.ac_org_y           = fread(fid, 1,  'int16',   0, 'b');
data.ac_width           = fread(fid, 1,  'int16',   0, 'b');
data.ac_height          = fread(fid, 1,  'int16',   0, 'b');
data.ac_n_buckets       = fread(fid, 1,  'int16',   0, 'b');
data.ac_range           = fread(fid, 1,  'int16',   0, 'b');
data.ac_n_bytes         = fread(fid, 1,  'int32',   0, 'b');
data.cn_org_x           = fread(fid, 1,  'int16',   0, 'b');
data.cn_org_y           = fread(fid, 1,  'int16',   0, 'b');
data.cn_width           = fread(fid, 1,  'int16',   0, 'b');
data.cn_height          = fread(fid, 1,  'int16',   0, 'b');
data.cn_n_bytes         = fread(fid, 1,  'int32',   0, 'b');
data.time_stamp         = fread(fid, 1,  'int32',   0, 'b');
data.comment            = fread(fid, 82, '*char'         )';
data.source             = fread(fid, 1,  'int16',   0, 'b');
data.intf_scale_factor  = fread(fid, 1,  'float32', 0, 'b');
data.wavelength_in      = fread(fid, 1,  'float32', 0, 'b');
data.num_aperture       = fread(fid, 1,  'float32', 0, 'b');
data.obliquity_factor   = fread(fid, 1,  'float32', 0, 'b');
data.magnification      = fread(fid, 1,  'float32', 0, 'b');
data.lateral_res        = fread(fid, 1,  'float32', 0, 'b');
data.acq_type           = fread(fid, 1,  'int16',   0, 'b');
data.intens_avg_cnt     = fread(fid, 1,  'int16',   0, 'b');
data.ramp_cal           = fread(fid, 1,  'int16',   0, 'b');
data.sfac_limit         = fread(fid, 1,  'int16',   0, 'b');
data.ramp_gain          = fread(fid, 1,  'int16',   0, 'b');
data.part_thickness     = fread(fid, 1,  'float32', 0, 'b');
data.sw_llc             = fread(fid, 1,  'int16',   0, 'b');
data.target_range       = fread(fid, 1,  'float32', 0, 'b');
data.radCrvMeasureSeq   = fread(fid, 1,  'int16',   0, 'l');
data.min_mod            = fread(fid, 1,  'int32',   0, 'b');
data.min_mod_count      = fread(fid, 1,  'int32',   0, 'b');
data.phase_res          = fread(fid, 1,  'int16',   0, 'b');
data.min_area           = fread(fid, 1,  'int32',   0, 'b');
data.discon_action      = fread(fid, 1,  'int16',   0, 'b');
data.discon_filter      = fread(fid, 1,  'int32',   0, 'b');
data.connect_order      = fread(fid, 1,  'int16',   0, 'b');
data.sign               = fread(fid, 1,  'int16',   0, 'b');
data.camera_width       = fread(fid, 1,  'int16',   0, 'b');
data.camera_height      = fread(fid, 1,  'int16',   0, 'b');
data.sys_type           = fread(fid, 1,  'int16',   0, 'b');
data.sys_board          = fread(fid, 1,  'int16',   0, 'b');
data.sys_serial         = fread(fid, 1,  'int16',   0, 'b');
data.inst_id            = fread(fid, 1,  'int16',   0, 'b');
data.obj_name           = fread(fid, 12, '*char'         )';
data.part_name          = fread(fid, 40, '*char'         )';
data.codev_type         = fread(fid, 1,  'int16',   0, 'b');
data.phase_avg_count    = fread(fid, 1,  'int16',   0, 'b');
data.sub_sys_err        = fread(fid, 1,  'int16',   0, 'b');
                          fread(fid, 16, 'int8'           ); % not used
data.part_ser_num       = fread(fid, 40, '*char'         )';
data.refractive_index   = fread(fid, 1,  'float32', 0, 'b');
data.rem_tilt_bias      = fread(fid, 1,  'int16',   0, 'b');
data.rem_fringes        = fread(fid, 1,  'int16',   0, 'b');
data.max_area           = fread(fid, 1,  'int32',   0, 'b');
data.setup_type         = fread(fid, 1,  'int16',   0, 'b');
data.wrapped            = fread(fid, 1,  'int16',   0, 'b');
data.pre_connect_filer  = fread(fid, 1,  'float32', 0, 'b');
data.wavelength_in_2    = fread(fid, 1,  'float32', 0, 'b');
data.wavelength_fold    = fread(fid, 1,  'int16',   0, 'b');
data.wavelength_in_1    = fread(fid, 1,  'float32', 0, 'b');
data.wavelength_in_3    = fread(fid, 1,  'float32', 0, 'b');
data.wavelength_in_4    = fread(fid, 1,  'float32', 0, 'b');
data.wavelen_select     = fread(fid, 8,  '*char'         )';
data.fda_res            = fread(fid, 1,  'int16',   0, 'b');
data.scan_descr         = fread(fid, 20, '*char'         )';
data.n_fiducials        = fread(fid, 1,  'int16',   0, 'b');
data.fiducials          = fread(fid, 14, 'float32', 0, 'b');
data.pixel_width        = fread(fid, 1,  'float32', 0, 'b');
data.pixel_height       = fread(fid, 1,  'float32', 0, 'b');
data.exit_pupil_diam    = fread(fid, 1,  'float32', 0, 'b');
data.light_level_pct    = fread(fid, 1,  'float32', 0, 'b');
data.coords_state       = fread(fid, 1,  'int32',   0, 'l');
data.coords.x_pos       = fread(fid, 1,  'float32', 0, 'l');
data.coords.y_pos       = fread(fid, 1,  'float32', 0, 'l');
data.coords.z_pos       = fread(fid, 1,  'float32', 0, 'l');
data.coords.x_rot       = fread(fid, 1,  'float32', 0, 'l');
data.coords.y_rot       = fread(fid, 1,  'float32', 0, 'l');
data.coords.z_rot       = fread(fid, 1,  'float32', 0, 'l');
data.coherence_mode     = fread(fid, 1,  'int16',   0, 'l');
data.surface_filter     = fread(fid, 1,  'int16',   0, 'l');
data.sys_err_file_name  = fread(fid, 28, '*char'         )';
data.zoom_descr         = fread(fid, 8,  '*char'         )';
data.alpha_part         = fread(fid, 1,  'float32', 0, 'l');
data.beta_part          = fread(fid, 1,  'float32', 0, 'l');
data.dist_part          = fread(fid, 1,  'float32', 0, 'l');
data.cam_split_loc_x    = fread(fid, 1,  'int16',   0, 'l');
data.cam_split_loc_y    = fread(fid, 1,  'int16',   0, 'l');
data.cam_split_trans_x  = fread(fid, 1,  'int16',   0, 'l');
data.cam_split_trans_y  = fread(fid, 1,  'int16',   0, 'l');
data.material_a         = fread(fid, 24, '*char'         )';
data.material_b         = fread(fid, 24, '*char'         )';
data.cam_split_unused   = fread(fid, 1,  'int16',   0, 'l');
                          fread(fid, 1,  'int16',   0, 'l'); % not used
data.dmi_ctr_x          = fread(fid, 1,  'float32', 0, 'l');
data.dmi_ctr_y          = fread(fid, 1,  'float32', 0, 'l');
data.sph_dist_corr      = fread(fid, 1,  'int16',   0, 'l');
                          fread(fid, 1,  'int16',   0, 'l'); % not used
data.sph_dist_part_na   = fread(fid, 1,  'float32', 0, 'l');
data.sph_dist_part_radius=fread(fid, 1,  'float32', 0, 'l');
data.sph_dist_cal_na    = fread(fid, 1,  'float32', 0, 'l');
data.sph_dist_cal_radius= fread(fid, 1,  'float32', 0, 'l');
data.surface_type       = fread(fid, 1,  'int16',   0, 'l');
data.ac_surface_type    = fread(fid, 1,  'int16',   0, 'l');
data.zPosition          = fread(fid, 1,  'float32', 0, 'l');
data.powerMultiplier    = fread(fid, 1,  'float32', 0, 'l');
data.focusMultiplier    = fread(fid, 1,  'float32', 0, 'l');
data.radCrvFocusCalFactor=fread(fid, 1,  'float32', 0, 'l');
data.radCrvPowerCalFactor=fread(fid, 1,  'float32', 0, 'l');
data.ftp_left_pos       = fread(fid, 1,  'float32', 0, 'l');
data.ftp_right_pos      = fread(fid, 1,  'float32', 0, 'l');
data.ftp_pitch_pos      = fread(fid, 1,  'float32', 0, 'l');
data.ftp_roll_pos       = fread(fid, 1,  'float32', 0, 'l');
data.min_mod_pct        = fread(fid, 1,  'float32', 0, 'l');
data.max_inten          = fread(fid, 1,  'int32',   0, 'l');
data.ring_of_fire       = fread(fid, 1,  'int16',   0, 'l');
                          fread(fid, 1,  'int8',    0, 'l'); % not used
data.rc_orientation     = fread(fid, 1,  '*char'         )';
data.rc_distance        = fread(fid, 1,  'float32', 0, 'l');
data.rc_angle           = fread(fid, 1,  'float32', 0, 'l');
data.rc_diameter        = fread(fid, 1,  'float32', 0, 'l');
data.rem_fringes_mode   = fread(fid, 1,  'int16',   0, 'b');
                          fread(fid, 1,  'int8'          ); % not used
data.ftpsi_phase_res    = fread(fid, 1,  'int8'           );
data.frames_acquired    = fread(fid, 1,  'int16',   0, 'l');
data.cavity_type        = fread(fid, 1,  'int16',   0, 'l');
data.cam_frame_rate     = fread(fid, 1,  'float32', 0, 'l');
data.tune_range         = fread(fid, 1,  'float32', 0, 'l');
data.cal_pix_loc_x      = fread(fid, 1,  'int16',   0, 'l');
data.cal_pix_loc_y      = fread(fid, 1,  'int16',   0, 'l');
data.n_tst_cal_pts      = fread(fid, 1,  'int16',   0, 'l');
data.n_ref_cal_pts      = fread(fid, 1,  'int16',   0, 'l');
data.tst_cal_pts        = fread(fid, 4,  'float32', 0, 'l');
data.ref_cal_pts        = fread(fid, 4,  'float32', 0, 'l');
data.tst_cal_pix_opd    = fread(fid, 1,  'float32', 0, 'l');
data.ref_cal_pix_opd    = fread(fid, 1,  'float32', 0, 'l');
data.flash_phase_cd_mask =fread(fid, 1,  'float32', 0, 'l');
data.flash_phase_alias_ma=fread(fid, 1,  'float32', 0, 'l');
data.flash_ogase_filter = fread(fid, 1,  'float32', 0, 'l');
data.scan_direction     = fread(fid, 1,  'int8'           );
                          fread(fid, 3,  'int8'           ); % not used
data.ftpsi_res_factor   = fread(fid, 1,  'int32',   0, 'l');
                          fread(fid, 16, 'int8'           ); % not used

% for MetroPro versions < 8.1.1 the header size is 834 bytes
% for higher versions the header size is 4096 bytes
fseek(fid, data.header_size, -1);
% if data.swinfo.vers.maj <= 8 & data.swinfo.vers.min <= 1 & data.swinfo.vers.bug < 1
%     fseek(fid, 834, -1);
% else                    
%     fseek(fid, 4096, -1);
% end


%
%%% INTENSITY DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
BAD = 65535;
ndata = data.ac_width * data.ac_height * data.ac_n_buckets;
data.IntensityData = [];

if ndata>0
    [tmp, count] = fread(fid, ndata, 'uint16', 0, 'b');
    if count<ndata
        fclose(fid);
        error('corrupted file in intensity data');
    end
    
    data.IntensityData = reshape(tmp, data.ac_width, data.ac_height, data.ac_n_buckets);
    idx = data.IntensityData >= BAD;
    data.IntensityData(idx) = NaN;
end

%
%%% PHASE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
BAD = 2147183640;
ndata = data.cn_width * data.cn_height;
data.PhaseData = [];

if ndata>0
    [tmp, count] = fread(fid, ndata, 'int32', 0, 'b');
    if count<ndata
        fclose(fid);
        error('corrupted file in phase data');
    end
    
    data.PhaseData = reshape(tmp, data.cn_width, data.cn_height);
    idx = data.PhaseData >= BAD;
    data.PhaseData(idx) = NaN;
end

fclose(fid);
