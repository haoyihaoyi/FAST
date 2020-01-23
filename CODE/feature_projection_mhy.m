function x = feature_projection_mhy(ori_feature,projection_matrix, rd_dim, w, h)
         x = reshape(ori_feature*projection_matrix, [w, h, rd_dim]);
end