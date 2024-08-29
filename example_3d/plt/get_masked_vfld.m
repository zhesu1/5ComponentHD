
function omega_vf_masked = get_masked_vfld(omega_vf, mask)

omega_vf_masked = zeros(size(omega_vf));
for i=1:3
    
    tem = omega_vf(:,:,:,i);
    tem(~mask) = 0;
    omega_vf_masked(:,:,:,i) = tem;

end
end