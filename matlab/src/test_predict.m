
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Change this

xp_dir = '/cvlabdata1/home/lucchi/xp/EM/striatum/autostep/training_files_striatum_0_1_5_f33_dist2_w9_c10_l0.001_t3_m0.0_r1e-05_i9_n0_o2_cs10000_cst0_sk1.0_rg0_s1/';
params = [xp_dir 'parameter_vector0/iteration_1440.txt'];
filename = '/cvlabdata1/home/biomed/EM/datasets/train_striatum/';
gt_filename = '/cvlabdata1/cvlab/datasets_carlos/striatum/train/mito_masks/';
config_file = [xp_dir '/training_files_striatum_0_1_5_f33_dist2.txt'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


close all;

addpath('../mex/bin');

labels = mex_predict(config_file, params, filename);

figure;

coeff_mul = floor(size(labels,1)/4);;
for i=1:4
    subplot(2,2,i);
    img_idx = i*coeff_mul;
    imshow(squeeze(labels(img_idx,:,:))*255);
    title(img_idx);
end
