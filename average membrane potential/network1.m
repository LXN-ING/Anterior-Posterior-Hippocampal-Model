function net = network1(n_ext)

% 初始化结构体
net = struct();

% 初始化 A 群体的连接
net.neigh_send_internal_A = zeros(500, 50);
net.neigh_receiv_internal_A = zeros(500, 50);
net.neigh_receiv_external_A = zeros(500, 20);
net.neigh_send_external_A = zeros(500, 20);
net.neigh_A2_from_P5 = zeros(500, n_ext);
net.neigh_A5_from_P5 = zeros(500, n_ext);
net.neigh_A2_from_P2 = zeros(500, n_ext);
net.neigh_A5_from_P2 = zeros(500, n_ext);

for j = 1:500
    for jj = 1:50
        net.neigh_send_internal_A(j, jj) = randperm(500, 1);
        if net.neigh_send_internal_A(j, jj) == j
            while net.neigh_send_internal_A(j, jj) == j
                net.neigh_send_internal_A(j, jj) = randperm(500, 1);
            end
        end
    end
end

for j = 1:500
    for jj = 1:50
        net.neigh_receiv_internal_A(j, jj) = randperm(500, 1);
        if net.neigh_receiv_internal_A(j, jj) == j
            while net.neigh_receiv_internal_A(j, jj) == j
                net.neigh_receiv_internal_A(j, jj) = randperm(500, 1);
            end
        end
    end
end

for j = 1:500
    for jj = 1:20
        net.neigh_receiv_external_A(j, jj) = randperm(400, 1);
        net.neigh_send_external_A(j, jj) = randperm(400, 1);
    end
end

for j = 1:500
    for jj = 1:n_ext
        net.neigh_A2_from_P5(j, jj) = randperm(400, 1);
        net.neigh_A5_from_P5(j, jj) = randperm(400, 1);
        net.neigh_A2_from_P2(j, jj) = randperm(400, 1);
        net.neigh_A5_from_P2(j, jj) = randperm(400, 1);
    end
end

% 初始化 P 群体的连接
net.neigh_send_internal_P = net.neigh_send_internal_A;
net.neigh_receiv_internal_P = net.neigh_receiv_internal_A;
net.neigh_receiv_external_P = net.neigh_receiv_external_A;
net.neigh_send_external_P = net.neigh_send_external_A;

net.neigh_P2_from_A5 = net.neigh_A2_from_P5;
net.neigh_P5_from_A5 = net.neigh_A5_from_P5;
net.neigh_P2_from_A2 = net.neigh_A2_from_P2;
net.neigh_P5_from_A2 = net.neigh_A5_from_P2;

end
