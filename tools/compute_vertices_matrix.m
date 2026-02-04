function IS_vertices_matrix = compute_vertices_matrix(IS)
% Given the intervalMatirx IS (CORA object) computes the vertices and
% stucks them in a 3D array

    [n,nm] = size(IS);
    IS_vertices = vertices(IS);
    nverts = length(IS_vertices);
    IS_vertices_matrix = zeros(n,nm,nverts);
    for i = 1:nverts
        IS_vertices_matrix(:,:,i) = IS_vertices{i};
    end
end