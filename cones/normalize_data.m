function [data, w] = normalize_data(data, K, scale, w)
[m,n] = size(data.A);

MIN_SCALE = 1e-3;
MAX_SCALE = 1e3;
minRowScale = MIN_SCALE * sqrt(n);
maxRowScale = MAX_SCALE * sqrt(n);
minColScale = MIN_SCALE * sqrt(m);
maxColScale = MAX_SCALE * sqrt(m);

D = ones(m,1);
E = ones(n,1);
NN = 1; % NN = 1, other choices bad
for j=1:NN
    %% D scale:
    Dt = twonorms(data.A(1:K.zero,:)')';
    idx = K.zero;
    Dt = [Dt;twonorms(data.A(idx+1:idx+K.nn,:)')'];
    idx = idx + K.nn;
    for i=1:length(K.soc)
        if (K.soc(i) > 0)
            nmA = mean(twonorms(data.A(idx+1:idx+K.soc(i),:)'));
            Dt = [Dt;nmA*ones(K.soc(i),1)];
            idx = idx + K.soc(i);
        end
    end
%     for i=1:length(K.sdc)
%         if (K.sdc(i) > 0)
%             nmA = mean(twonorms(data.A(idx+1:idx+getSdConeSize(K.sdc(i)),:)'));
%             Dt = [Dt;nmA*ones(getSdConeSize(K.sdc(i)),1)];
%             idx = idx + getSdConeSize(K.sdc(i));
%         end
%     end
%     for i=1:K.expp
%         nmA = mean(twonorms(data.A(idx+1:idx+3,:)'));
%         Dt = [Dt;nmA*ones(3,1)];
%         idx = idx + 3;
%     end
%     for i=1:K.expd
%         nmA = mean(twonorms(data.A(idx+1:idx+3,:)'));
%         Dt = [Dt;nmA*ones(3,1)];
%         idx = idx + 3;
%     end
%     for i=1:length(K.pow)
%         nmA = mean(twonorms(data.A(idx+1:idx+3,:)'));
%         Dt = [Dt;nmA*ones(3,1)];
%         idx = idx + 3;
%     end
    
    Dt(Dt < minRowScale) = 1;
    Dt(Dt > maxRowScale) = maxRowScale;
    data.A = sparse(diag(1./Dt))*data.A;
    
    %% E Scale
    Et = twonorms(data.A)';
    Et(Et < minColScale) = 1;
    Et(Et > maxColScale) = maxColScale;
    data.A = data.A*sparse(diag(1./Et));
    
    %%
    D = D.*Dt;
    E = E.*Et;
end

nmrowA = mean(twonorms(data.A'));
nmcolA = mean(twonorms(data.A));

data.A = data.A*scale;

data.b = data.b./D;
sc_b = nmcolA/ max(norm(data.b), MIN_SCALE);
data.b = data.b * sc_b * scale;

data.c = data.c./E;
sc_c = nmrowA/max(norm(data.c), MIN_SCALE);
data.c = data.c * sc_c * scale;

w.D = D;
w.E = E;
w.sc_b = sc_b;
w.sc_c = sc_c;

    function twoNorms = twonorms(A)
        twoNorms = sqrt(sum(A.^2,1));
    end

end