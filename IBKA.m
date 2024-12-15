
%%  Black-winged Kite Algorithm

function [Best_Fitness_BKA,Best_Pos_BKA,Convergence_curve]=IBKA(SearchAgents_no,Max_iter,lb,ub,dimension,fobj)
%% ----------------Initialize the locations of Blue Sheep------------------%
p=0.9;
r=rand; 
lammda = 0.4;
Bernoulli=rand(SearchAgents_no,dimension); 
for i=1:SearchAgents_no
    for j=2:dimension
        if Bernoulli(i,j-1) <  1-lammda
            Bernoulli(i,j)= Bernoulli(i,j-1)/(1-lammda);
        else
            Bernoulli(i,j)= (Bernoulli(i,j-1)-1+lammda)/lammda;
        end
    end
end
result = Bernoulli;

XPos=repmat(lb,SearchAgents_no,1)+result.* repmat((ub-lb),SearchAgents_no,1);


for i =1:SearchAgents_no
    XFit(i)=fobj(XPos(i,:));
end
Convergence_curve=zeros(1,Max_iter);
%% -------------------Start iteration------------------------------------%
for t=1:Max_iter
    p=rand; 
    [~,sorted_indexes]=sort(XFit);
    XLeader_Pos=XPos(sorted_indexes(1),:);
    XLeader_Fit = XFit(sorted_indexes(1));
   
%% -------------------Attacking behavior-------------------%

        for i=1:SearchAgents_no
        if rand>0.5
        n=0.05*exp(-2*(t/Max_iter)^2);
        if p<r
            XPosNew(i,:)=XPos(i,:)+n.*(1+sin(r))*XPos(i,:);
        else
            XPosNew(i,:)= XPos(i,:).*(n*(2*rand(1,dimension)-1)+1);
        end
        else
         R=0.02*(1-t/Max_iter);% Eq.(6)
        XPosNew(i,:)= XPos(i,:)+ (-R+2*R*rand(1,dimension)).*XPos(i,:);% Eq.(7)
        end        

        XPosNew(i,:) = max(XPosNew(i,:),lb);
        XPosNew(i,:) = min(XPosNew(i,:),ub);%%Boundary checking
%% ------------ Select the optimal fitness value--------------%
        
        XFit_New(i)=fobj(XPosNew(i,:));
        if(XFit_New(i)<XFit(i))
            XPos(i,:) = XPosNew(i,:);
            XFit(i) = XFit_New(i);
        end
%% -------------------Migration behavior-------------------%
        m=2*sin(r+pi/2);
        s = randi([1,30],1);
        r_XFitness=XFit(s);
        ori_value = rand(1,dimension);
        cauchy_value = tan((ori_value-0.5)*pi);
         R=0.02*(1-t/Max_iter);% Eq.(6)
       % Randomly select 2 to 5 red-billed blue magpies
         p = randi([2, 5]);
         selected_index_p = randperm(SearchAgents_no, p);
         Xp = XPos(selected_index_p, :);
         Xpmean = mean(Xp);

        % Randomly select 10 to N red-billed blue magpies
         q = randi([10, SearchAgents_no]);
         selected_index_q = randperm(SearchAgents_no, q);
         Xq = XPos(selected_index_q, :);
         Xqmean = mean(Xq);
        if XFit(i)< r_XFitness
            XPosNew(i,:)=XPos(i,:)+cauchy_value(:,dimension).* (XPos(i,:)-Xpmean);  
        else
            XPosNew(i,:)=XPos(i,:)+cauchy_value(:,dimension).* (Xqmean-m.*XPos(i,:));  
        end
        XPosNew(i,:) = max(XPosNew(i,:),lb);XPosNew(i,:) = min(XPosNew(i,:),ub); %%Boundary checking
        %% --------------  Select the optimal fitness value---------%
        XFit_New(i)=fobj(XPosNew(i,:));
        if(XFit_New(i)<XFit(i))
            XPos(i,:) = XPosNew(i,:);
            XFit(i) = XFit_New(i);
        end
%% Secretary Bird's escape strategy     
            Epsilon = unifrnd(0, 1);
            F=0.5;
            rand_leader_index1 = floor(SearchAgents_no * rand() + 1);
            rand_leader_index2 = floor(SearchAgents_no * rand() + 1);
            X_rand1 = XPos(rand_leader_index1, :);
            X_rand2 = XPos(rand_leader_index2, :);
            step2 = X_rand1 - X_rand2;
            if rand<0.5
                XPosNew(i,:) = XPos(i,:) + Epsilon .* step2; 
            else
                XPosNew(i,:) = XPos(i,:) + F .* Levy(dimension) .* step2; 
            end    
        XPosNew(i,:) = max(XPosNew(i,:),lb);XPosNew(i,:) = min(XPosNew(i,:),ub); %%Boundary checking
        %% --------------  Select the optimal fitness value---------%
        XFit_New(i)=fobj(XPosNew(i,:));
        if(XFit_New(i)<XFit(i))
            XPos(i,:) = XPosNew(i,:);
            XFit(i) = XFit_New(i);
        end

    end
    %% -------Update the optimal Black-winged Kite----------%
    if(XFit<XLeader_Fit)
        Best_Fitness_BKA=XFit(i);
        Best_Pos_BKA=XPos(i,:);
    else
        Best_Fitness_BKA=XLeader_Fit;
        Best_Pos_BKA=XLeader_Pos;
    end
    Convergence_curve(t)=Best_Fitness_BKA;
end
end
%% Levy flight function
function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end