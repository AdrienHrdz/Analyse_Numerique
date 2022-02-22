clear variables;
close all;
clc;

exo = 3;

switch exo
    case 1
        %% Exercice 1 
        A = [2,2 ; 1,3];
        B = inv(A);
        norme = 1;
        switch norme
            case 1
               K = norm(A,1)*norm(B,1);
               K_cond = cond(A,1);
            case 2
                K = norm(A,2)*norm(B,2);
                K_cond = cond(A,2);
            case 3
                K = norm(A,"inf")*norm(B,"inf");
                K_cond = cond(A,"inf");
        end
    case 2
        %% Exercice 2 
        A = [8 10 10 10;
            1 5 6 1;
            7 10 4 7;
            7 8 1 7];
        b = [38 13 28 23]';

        x = A\b;
        
        db = [0.1 -0.1 0.1 -0.1]';

        x1 = A\(b + db);

        dx = x1 - x;
        qb = norm(db)/norm(b);
        qx = norm(dx)/norm(x);
        
        S = svd(A);
        K_2 = S(1)/S(length(S));
        K_2_Cond = cond(A);
        if( qx <= K_2*qb )
            disp("formule 4 du cours vérifiée")
        end
        
        epsi = 5e-4;
        dA = [epsi 0 0 epsi;
            0 epsi 0 epsi;
            0 epsi epsi 0;
            epsi 0 0 epsi];

        x2 = (A + dA)\b;
        dx2 = x2 - x;
        qx2 = norm(dx2)/norm(dx2 + x2);
        qA = norm(dA)/norm(A);
        if ( qx2 <= K_2*qA )
            disp('formule 5 du cours vérifiée')
        end

    case 3
        list_produit = zeros(1,9);
        for n = 4:12
            A = zeros(n,n);
            for i = 1:n
                for j = 1:n
                    A(i,j) = 1/(i + j -1);
                end
            end
            K = cond(A);
            produit = det(A)*det(inv(A));
            list_produit(n-3) = produit;
        end
        % on remarque que le produit théorique n'est plus rigoureusement
        % égale à 1 lorsque la matrice se rapproche de la matrice identité

        A = hilb(5);
        b = [137/60 29/20 153/140 743/840 1879/2520]';

        [U, S, V] = svd(A);
        x = A\b;
        db = 1/1000*[0 1 0 1 0]';
        x1 = A\(b + db);
        K = cond(A);
        % on trouve un résultat totalement erroné par rapport à la réponse attendue
        % cela s'explique par le fait que la matrice A est mal conditionnée
        % en effet K(A) est environ égale à 500 000.

        lambda = 0.001;
        x_tik_formule8 = inv(A'*A + lambda*eye(size(A)))*A'*b;

        r = rank(A);
        x_tik_formule9 = zeros(size(x));
        for k = 1:r
            x_tik_formule9 = x_tik_formule9 + (S(k,k)/(S(k,k)^2+lambda)) * V(:,k) * U(:,k)' * b;
        end
        % avec les formules 8 et 9 du cours, on retrouve le même résultat
        % pour x. En effet ces formules sont identiques, seule la façon de
        % les coder diffère. 
        % De plus, on peut voir que le résultat obtenu à l'aide de la
        % régularisation est nettement plus proche du résultat théorique.
end