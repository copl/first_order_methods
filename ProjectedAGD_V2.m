function [ x, y, s ] = ProjectedAGD_V2(c, A, b, it, eps )
    % solves Ax = b, x >= 0
    % by solving (x+)^2
    %newton_matrix = [2*eye(size(A,2)) A';
    %A zeros(size(A,1))];
    %value_indicies = [];
    %values = [];
    %store_index = 1;
    %alpha = 1.1;
    
    %[L,U] = lu(newton_matrix);
    t = 1.0/4.0;
    tic
        [R_primal,p_primal,S_primal] = chol(sparse(A*A'+ 10^(-8.0)*eye(size(A,1)))); % + 10^-7*eye(size(A,1))));
        assert(p_primal == 0) % is matrix positive definite??? 
        
        [R_dual,p_dual,S_dual] = chol(sparse(A'*A + eye(size(A,2))));
        assert(p_dual == 0) % is matrix positive definite??? 
    toc
    
    [m,n] = size(A);
    x = zeros(n,1);
    y = zeros(m,1);
    s = zeros(n,1);
    xOLD = x;
    yOLD = y;
    sOLD = s;
    
    tic
        duality_gap_scale = max(max(c)*norm(c),max(b)*norm(b));
        for k = 0:it
            omega = k/(k + 3.0);
            xEval = x + omega*(x - xOLD);
            yEval = y + omega*(y - yOLD);
            sEval = s + omega*(s - sOLD);
            
            %  
            xGrad = 2*c*(c'*xEval - b'*yEval)/duality_gap_scale + 2*min(xEval,0);
            xSuggest = xEval - t*xGrad;
            
            yGrad = -2*b*(c'*xEval - b'*yEval)/duality_gap_scale;
            ySuggest = yEval - t*yGrad;
            
            sGrad = 4*min(sEval,0);
            sSuggest = sEval - t*sGrad;
                       
            % project onto primal
            primalRhs = 2*(A*xSuggest - b);
            primalPi = (S_primal * (R_primal \ (R_primal' \ (S_primal' * primalRhs)))); % linsolve
            x = xSuggest - 0.5*(A'*primalPi);
            
            % project onto dual            
            dualRhs = 2*((A'*ySuggest) - sSuggest - c);
            dualPi = (S_dual * (R_dual \ (R_dual' \ (S_dual' * dualRhs))));
            y = ySuggest - 0.5*A*dualPi;
            s = sSuggest + 0.5*dualPi;
            
            %norm(min(x,0))
            %norm(min(s,0))
            
            xOLD = x;
            yOLD = y;
            sOLD = s;
            % if k >= store_index 
            %     values = [values, x];
            %     value_indicies = [value_indicies, k];
            %     store_index = store_index*alpha;
            % end

            if any(isnan(x))
                disp('error NaN in x')
                break
            end

            %if obj < eps
            %    disp('eps feasible solution found')
            %    break
            %end

            if k == it
                %norm(grad,2)
                disp('iteration limit reached')
            end
        end
        
        %values = [values, x];
        %value_indicies = [value_indicies, k];
        %store_index = store_index*alpha;
    toc
end

