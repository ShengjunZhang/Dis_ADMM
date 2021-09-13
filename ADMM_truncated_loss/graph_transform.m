function [Adj, degree, num_of_edge,A,B,D,Lm,edge_index, eig_Lm,min_eig_Lm,WW,LN,L_hat,eig_L_hat,min_eig_L_hat] = graph_transform(Adj, n_agents, d)

    degree=Adj*ones(n_agents,1);

    num_of_edge = sum(sum(Adj))/2;
    A = zeros(num_of_edge, n_agents);
    edge_index = 1;
    for ii = 1 : n_agents
        for jj = 1 : n_agents
            if (jj>ii)
                if Adj(ii,jj) == 1
                    A(edge_index,jj) = 1;
                    A(edge_index,ii) = -1;
                    edge_index = edge_index + 1;
                end
            end
        end
    end
    A = kron(A,eye(d));
    B = abs(A);
    Lm = A' * A; % signed LAplacian
    D = (B' * B + A' * A)/2;


    eig_Lm = eig(Lm);
    for ii = 1 : length(eig_Lm)
        if (eig_Lm(ii)>=1e-10)
            min_eig_Lm = eig_Lm(ii);
            break;
        end
    end

    WW = zeros(size(Adj));
    for ii = 1:size(Adj,1)
        for jj = 1:size(Adj,1)
            if Adj(ii, jj) == 1
                WW(ii,ii) = WW(ii,ii)+1/sqrt(degree(ii)*degree(jj));
            end
        end
    end
    WW = kron(WW, eye(d));

    LN = D^(-1/2)*(A'*A)*D^(-1/2);
    L_hat = LN - diag(diag(LN)) + diag(diag(WW));
    eig_L_hat = eig(L_hat);

    for ii = 1 : length(eig_L_hat)
        if (eig_L_hat(ii)>=1e-10)
            min_eig_L_hat = eig_L_hat(ii);
            break;
        end

    end
end