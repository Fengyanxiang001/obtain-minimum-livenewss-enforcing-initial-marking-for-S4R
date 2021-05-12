%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% author: Feng Yanxiang
% 2021/02/12
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clear
clc
MaxPro = 10;                         %Maximum number of process types
MaxStag = 8;                         %Maximum number of stages for each parocess type
MaxRes = 8;                          %Maximum number of resource types
MaxResUse = 8;                       %Maximum number of resources used by a stage
MaxResNu = 6;                        %Maximum number of resources types used by a stage
IM0Marking = 6;                      %Maximum number of token in each idel place initially

ProcesNu=randperm(MaxPro,1);        %the number of process types
ResourceNu=randperm(MaxRes,1);      %the number of resource types
StageNu=zeros(ProcesNu,2);          %the number of stages for each process
RouteNu=zeros(ProcesNu,1);          %the number of routes 
CommStag=zeros(ProcesNu,2);         %the ccommon stages of each process type   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%establish the S4R model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:1:ProcesNu
    RouteNu(i)=randperm(2,1);          
    if(RouteNu(i)>1)                  
        for j = 1:1: RouteNu(i)
            StageNu(i,j)=randperm(ceil(MaxStag/2),1) + ceil(MaxStag/2);     
        end      
        CommStag(i,1)=randperm(ceil(min(StageNu(i,1),StageNu(i,2))/4),1);    
        CommStag(i,2)=randperm(ceil(min(StageNu(i,1),StageNu(i,2))/4),1);
    else                                 
        StageNu(i,1)=randperm(MaxStag,1);
    end
end
NumberStage = 0;                       %the sum of all stages initially set as 0
for i=1:1:ProcesNu
    for j = 1:1:2
        NumberStage = StageNu(i,j)+NumberStage;
    end
end

for i=1:1:NumberStage
    while(1)
        ResUseNu(i,:) = floor(unifrnd(0,MaxResUse,1,ResourceNu));
        if(sum(ResUseNu(1,:)))>0
            break
        end
    end
end

TranNu=0;                                         %the number of transitions
ActPlaNu=0;                                       %the number of places
ZTraNu=zeros(ProcesNu,1);                         %ZTraNu[i]is the number of transitions in the i-th process
ZPlacNu=zeros(ProcesNu,1);                        %ZPlacNu[i] is the number of places in the i-th process 

for i=1:1:ProcesNu
    ss=StageNu(i,2)-CommStag(i,1)-CommStag(i,2);   
        
    for k=1:1:StageNu(i,1)+1                         
        N(k,k) = -1; 
    end
    for k=1:1:StageNu(i,1)                         
        N(k,k+1) =1;
    end
    N(StageNu(i,1)+1,1)=1;
    if(ss>0)      
        N(StageNu(i,1)+2,CommStag(i,1)+1)=-1;
        for k =1:1:ss
            N(StageNu(i,1)+1+k,StageNu(i,1)+1+k) = 1;
        end
        for h =1:1:ss
            N(StageNu(i,1)+2+h,StageNu(i,1)+1+h) = -1;
        end
        N(StageNu(i,1)+2+ss,StageNu(i,1)-CommStag(i,2)+2)=1;                   
        for (g = 1:1:StageNu(i,1)+2+ss)
            for(d = 1:1:StageNu(i,1)+1+ss)
                NN(g+TranNu,d+ActPlaNu)=N(g,d);
            end
        end
        ZTraNu(i)=StageNu(i,1)+2+ss;
        ZPlacNu(i)=StageNu(i,1)+1+ss;
        TranNu = TranNu + StageNu(i,1)+2+ss;
        ActPlaNu = ActPlaNu + StageNu(i,1)+1+ss;
    else                                           
        for (g = 1:1:StageNu(i,1)+1)
            for(d = 1:1:StageNu(i,1)+1)
                NN(g+TranNu,d+ActPlaNu)=N(g,d);
            end
        end
        ZTraNu(i)=StageNu(i,1)+1;
        ZPlacNu(i)=StageNu(i,1)+1;
        TranNu = TranNu + StageNu(i,1)+1;
        ActPlaNu = ActPlaNu + StageNu(i,1)+1;
    end
end


PreTranSet=zeros(ActPlaNu,2);   %the set of pre-transitions of each place
ProTranSet=zeros(ActPlaNu,2);   %the set of pro-transitions of each place
for i=1:1:ActPlaNu    
    a=1;
    b=1;
    for j=1:1:TranNu
        if NN(j,i)==-1 && (PreTranSet(i,a)~= j)
            PreTranSet(i,a)=j;
            a=a+1;
        end
        if a==3 | b==3
            break;
        end
        if NN(j,i)==1 && ProTranSet(i,b)~= j
            ProTranSet(i,b)=j;
            b=b+1;
        end
       
    end
end

NNR=zeros(TranNu,ActPlaNu+ResourceNu);  %the incidence matrix of S4R
for i=1:1:TranNu                       
    for j=1:1:ActPlaNu
        NNR(i,j)=NN(i,j);
    end
end

while(1)
    hh=0;
    ActPlaUseResNu = zeros(ActPlaNu,ResourceNu);
    for i=1:1:ProcesNu
    ss=StageNu(i,2)-CommStag(i,1)-CommStag(i,2);
    for j =2:1:ZPlacNu(i)
        while(1)              
            ResUseN = floor(unifrnd(0,MaxResUse,1,ResourceNu));
            if(sum(ResUseN)>0)
                break
            end
        end                  
       for k =1:1:ResourceNu  
           ActPlaUseResNu(j+hh,k) =  ResUseN(k);
           if ProTranSet(j+hh,2) == 0
               new_C=ProTranSet(j+hh,1);
               NNR(new_C,ActPlaNu+k)=NNR(new_C,ActPlaNu+k)-ResUseN(k);
           end
           if ProTranSet(j+hh,2) > 0 & ResUseN(k)>0
               aa1=ProTranSet(j+hh,1);
               aa2=ProTranSet(j+hh,2);
               NNR(aa1,ActPlaNu+k)=NNR(aa1,ActPlaNu+k)-ResUseN(k);
               NNR(aa2,ActPlaNu+k)=NNR(aa2,ActPlaNu+k)-ResUseN(k);               
           end
           if PreTranSet(j+hh,2) == 0
               bb=PreTranSet(j+hh,1);
               NNR(bb,ActPlaNu+k)=NNR(bb,ActPlaNu+k)+ResUseN(k);
           end
           if PreTranSet(j+hh,2) > 0 && ResUseN(k)>0
               bb1=PreTranSet(j+hh,1);
               bb2=PreTranSet(j+hh,2);
               NNR(bb1,ActPlaNu+k)=NNR(bb1,ActPlaNu+k)+ResUseN(k);
               NNR(bb2,ActPlaNu+k)=NNR(bb2,ActPlaNu+k)+ResUseN(k); 
           end
       end
    end
    hh=hh+ZPlacNu(i);
    end
    kk = 1;
    for s = 1:1:ResourceNu
        if sum(ActPlaUseResNu(:,s))==0
            kk=0;
            break;
        end
    end 
    if kk == 1
        break;
    end
end

M0=floor(unifrnd(2,IM0Marking,1,ProcesNu));    %the initial state of each idel place


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%compute all PACs 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M_NNR = zeros(TranNu,ResourceNu);               %the transition-resource subnet
for i = 1:1:TranNu
    for j = 1:1:ResourceNu
        M_NNR(i,j) = NNR(i,j+ActPlaNu);
    end
end

ANN = zeros(TranNu,ResourceNu);                %the virtual net
for i = 1:1:ActPlaNu
    a1 = PreTranSet(i,1);
    a2 = PreTranSet(i,2);
    for j = 1:1:ResourceNu
        if ActPlaUseResNu(i,j) >0 && a1 >0 && NNR(a1,ActPlaNu+j)~=1
            ANN(a1,j) = 1;
        end
        if ActPlaUseResNu(i,j) >0 && a2 >0 && NNR(a2,ActPlaNu+j)~=1
            ANN(a2,j) = 1;
        end
    end
end

BigN = zeros(TranNu+ResourceNu,TranNu+ResourceNu); 
for i =1:1:ResourceNu
    for j=1:1:TranNu
        if M_NNR(j,i)<0
            BigN(TranNu+i,j)=1;
        end
    end
end
for i =1:1:TranNu
    for j=1:1:ResourceNu
        if M_NNR(i,j)>0 || ANN(i,j)>0
            BigN(i,j+TranNu)=1;
        end
    end
end

%finding all elementary circuits 
ECircuitNu=0;   %the number of elementary circuits
ECircuit=[];    %the set of elementary circuits
for kkk = 1:1:ResourceNu   
    cE=TranNu+kkk;     
    V1=[];
    V1(1)=cE;          
    V2(1)=0;
    Flag(1)=0;
    current=1;
    PP=0;
    FlagB=0;
    FlagC=0;
    while(1)
        if FlagB ==1
            break;
        end
        if FlagC==0
           V=[];
           a=0;
           for j = 1:1:TranNu+ResourceNu 
               if BigN(cE,j)>0 
                   a=a+1;
                   V(a)=j;   
               end
           end
        else
            V=[];
        end    
        if length(V)==0 | FlagC==1   
            FlagA=0; 
            for i =1:1:current          
                if Flag(current-i+1)>0
                    current=current-i+1;
                    FlagA=1; 
                    for g=1:1:length(V2(current,:))
                        if ismember(V2(current,g),V1)>0
                            V2(current,g)=0;
                            Flag(current)=Flag(current)-1; 
                        end
                        if V2(current,g)>0 
                            cE=V2(current,g);
                            V2(current,g)=0;
                            Flag(current)=Flag(current)-1; 
                            break;
                        end
                    end  
                    for g=current+1:1:length(V1) 
                        V1(g)=0;
                    end                  
                    V1(current)=cE;                
                    break;
                end
            end
            if FlagA==0
                FlagB=1;
                break;
            end
            if FlagC==1
                FlagC=0;
            end
        end    
        if length(V)>0         
            FlagC = 0; 
            if ismember(V(1),V1)>0
                FlagC = 1;
                if V(1) == TranNu+kkk  
                    ECircuitNu=ECircuitNu+1;
                    for s=1:1:length(V1)
                        if V1(s)>0
                            ECircuit(ECircuitNu,s) = V1(s);
                        end
                    end
                   
                end
            end
            if FlagC ==0
                current=current+1;
                for j =1:1:length(V)-1
                    V2(current,j)=V(j+1);
                end
                Flag(current)=length(V)-1;
                V1(current)=V(1);
                cE=V1(current);
            end
        end
    end
end

%deleting all repeat elementary circuits
for i=1:1:ECircuitNu
    ECircuit(i,:)=sort(ECircuit(i,:),'descend');  %将所有回路序列按照从大到小重新排序
end
EC1=unique(ECircuit,'rows');
EC1_Num = size(EC1,1);

%compute all A-circuits
PEC=EC1;
for i=1:1:EC1_Num    
    j=1;
    while(j<=size(PEC,1))
        %fprintf('i=%d,j=%d\n',i,j);
        A1=EC1(i,:);
        A2=PEC(j,:);
        ss=[];
        ss=intersect(A1,A2);
        if length(find(ss>TranNu))>0 
            new_C=[];
            new_C=union(A1, A2);        
            new_C=sort(new_C,'descend');   
            ak=size(PEC,1);            
            F=0;
            for m=1:1:ak
                E=0;
                skk = min(length(new_C),size(PEC,2));
                for n=1:1:skk  
                    if PEC(m,n) ~= new_C(n)
                        E=1;
                        break;
                    end    
                end
                if E==0
                    F=1;
                    break;
                end
            end
            if F==0
                for t=1:1:length(new_C)
                    PEC(ak+1,t)=new_C(t);
                end
            end      
        end 
        j=j+1;
    end
end


%find all split places
splitTran=[];
aa=find(PreTranSet(:,2)>0);
for i =1:1:length(aa)
    splitTran(i,:) = PreTranSet(aa(i),:);
end

% find all perfect A-circuit,i.e., PAC
PC=[];
a=0;
for i=1:1:size(PEC,1)
    if size(splitTran,1) >0
        Gap=1;
        for j=1:1:size(splitTran,1)
            a1=ismember(splitTran(j,1),PEC(i,:));
            a2=ismember(splitTran(j,2),PEC(i,:));
            if a1+a2==1 
                Gap=0;
                break;
            end
        end
        if Gap==1
            a=a+1;
            PC(a,:)=PEC(i,:);
        end
    else
        PC=PEC;
    end
end


% comopute the transition set and resource set of each PAC
PC_R=[];
PC_T=[];
PC_P=[];
for i = 1:1:size(PC,1)
    aa=[];
    bb=[];
    aa=find(PC(i,:)>TranNu);
    bb=find(PC(i,:)<=TranNu & PC(i,:)>0) ;  
    for j=1:1:length(aa)
        PC_R(i,j)=PC(i,aa(j))-TranNu;
    end           
    for j=1:1:length(bb)
        PC_T(i,j)=PC(i,bb(j));    
        kk=PC_T(i,j);
        for k=1:1:ActPlaNu       
            if PreTranSet(k,1) == kk
                PC_P(i,j) = k;
                break
            end
            if PreTranSet(k,2) == kk
                PC_P(i,j) = k;
                break
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% compute all variables needed for obtaining a minimum LIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA=zeros(size(PC,1),ResourceNu);
MinW=zeros(size(PC,1),ResourceNu);
LL=zeros(size(PC,1),ResourceNu);
for i=1:1:size(PC,1)
    %fprintf('\n第%d条环\n',i);M
    for j=1:1:length(find(PC_R(i,:)>0))
        R=PC_R(i,j); 
        %fprintf('第%d个资源\n ',R);
        AA(i,R)=0;
        for k=1:1:ActPlaNu
            if ismember(k,PC_P(i,:))>0
                AA(i,R)=AA(i,R) + ActPlaUseResNu(k,R);
            end
        end                  
        be=1;
        LL(i,R)=0;
        for m=1:1:ProcesNu
            cc=[];
            edn=be+ZPlacNu(m)-1;
            cc=ActPlaUseResNu(be:edn,R);
            for p=be:1:edn           
                if ismember(p,PC_P(i,:))==0 && cc(p-be+1)>0
                    cc(p-be+1)=0;
                end
            end
            cc;
            be=edn+1;
            %fprintf('max(cc)=%d,min(abs(dd))=%d\n ',max(cc),min(abs(dd)));
            LL(i,R) = LL(i,R) + max(cc)*M0(m);
        end
        dd=find(NNR(:,ActPlaNu+R)<0);
        for p=1:1:length(dd)
            if ismember(dd(p),PC_T(i,:))==0 
                dd(p)=0;
            end
        end
        dd;  
        dd(find(dd==0))=[];
        dd1=[];
        for s=1:1:length(dd)
            dd1(s)=-NNR(dd(s),ActPlaNu+R);
        end
        if length(dd1)>0
            MinW(i,R)=min(dd1);
        else
            MinW(i,R)=0;
        end
        LL(i,R) = LL(i,R) + MinW(i,R);  %MinW(i,R)表示对应的权值的最小值
    end
end
AA;
LL;


BB=zeros(ResourceNu,1,size(PC,1));
BB_Num=[];
for i=1:1:size(PC,1)
    %fprintf('\n第%d条环\n',i);
    for j=1:1:length(find(PC_R(i,:)>0))
        R=PC_R(i,j); 
        PC_PR=[];
        PC_PR=PC_P(i,:);
        for k=1:1:length(PC_PR)
            if PC_PR(k)==0
                break;
            end
            if PC_PR(k)>0 
                currentPla= PC_PR(k);
                if ActPlaUseResNu(currentPla,R)==0 & PC>0
                    PC_PR(k)=0;
                end
            end
        end
        PC_PR(find(PC_PR==0))=[]; 
        a=0;
        for k=LL(i,R)-AA(i,R)-1:-1:1
            af=length(PC_PR);
            f=[];
            A1=zeros(af);
            B1=[];
            C1=zeros(af,1);
            D1=[];
            inc=[];
            for s=1:1:length(PC_PR)
                inc(s)=s;
                currentPla=PC_PR(s);
                f(s)=-ActPlaUseResNu(currentPla,R);
                D1(s)=-k/f(s);
                B1=[k];
                A1=-f;        
            end
            D1=D1';
            [x,fval]=intlinprog(f,inc,A1,B1,[],[],C1,D1);
            if(k+fval>=MinW(i,R))
                a=a+1;
                BB(R,a,i)=k;
            end
        end  
        BB_Num(i,R) = a;
    end
end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% output the results in the form of ".txt". The first file named by
% "Data_S4R" collects all information the randomly generated S4R model.
% The other "Data needed for computing Minimum LIM" collects
% all and variables needed for computing the corresponding LIM. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fis = fopen('Data_S4R.txt','wt');
fprintf(fis,'the number of process types: %d\n',ProcesNu);
fprintf(fis,'the umber of resources: %d\n',ResourceNu);
fprintf(fis,'the umber of flexible routes:\n');
for i=1:1:ProcesNu
        fprintf(fis,'%d ',RouteNu(i));
end
fprintf(fis,'\n the number of stages or operation for each process type:\n');
for i=1:1:ProcesNu
    for j=1:1:2
        fprintf(fis,'%d ',StageNu(i,j));
    end
    fprintf(fis,'\n');
end
fprintf(fis,'the number of shared operations for each process type:\n');
for i=1:1:ProcesNu
    for j=1:1:2
        fprintf(fis,'%d ',CommStag(i,j));
    end
    fprintf(fis,'\n');
end

fprintf(fis,'the initial state for each idle place is:\n');
for i=1:1:ProcesNu
    fprintf(fis,'%d ',M0(i));
end
fprintf(fis,'\n');

fprintf(fis,'the number of transitions of PN model: %d\n',TranNu);
fprintf(fis,'the number of places of PN model: %d\n',ActPlaNu+ResourceNu);
fprintf(fis,'the resource useage:\n');
for i=1:1:ActPlaNu
    for j=1:1:ResourceNu
        fprintf(fis,'%d ',ActPlaUseResNu(i,j));
    end
    fprintf(fis,'\n');
end

fprintf(fis,'the incidence matrix of PN model:\n');
for i=1:1:TranNu
    for j=1:1:ActPlaNu+ResourceNu
        fprintf(fis,'%d ',NNR(i,j));
    end
    fprintf(fis,'\n');
end

fprintf(fis,'the number of elementary circuits: %d\n',EC1_Num);
if EC1_Num >0 
    fprintf(fis,'all elementary circuits：\n');
    for i=1:1:size(EC1,1)
        for j=1:1:size(EC1,2)
            if EC1(i,j) >0
            fprintf(fis,'%d ',EC1(i,j));
            else
                break;
            end
        end
        fprintf(fis,'\n');
    end
end

fprintf(fis,'the number of PAC: %d\n',size(PC,1));
if size(PC,1) >0 
    fprintf(fis,'all PACs：\n');
    for i=1:1:size(PC,1)
        for j=1:1:size(PC,2)
            if PC(i,j)>0
            fprintf(fis,'%d ',PC(i,j));
            else
                break
            end
        end
        fprintf(fis,'\n');
    end
end




fclose(fis);


fid = fopen('Data needed for computing Minimum LIM.txt','wt');
fprintf(fid,'%d\n',ResourceNu);
fprintf(fid,'%d\n',size(BB(:,:,1),2));
fprintf(fid,'%d\n',size(PC,1));
fprintf(fid,'100000\n');
fprintf(fid,'10000\n');


fprintf(fid,'[');
for i= 1:1:size(AA,1)
    fprintf(fid,'[');
    for j = 1:1:ResourceNu
        if j < ResourceNu
            fprintf(fid,'%d,',AA(i,j));
        else
            fprintf(fid,'%d',AA(i,j));
        end
    end
    if i < length(AA)
        fprintf(fid,'],');
    else
        fprintf(fid,']');
    end
end
fprintf(fid,']');



fprintf(fid,'\n[');
for i= 1:1:length(LL)
    fprintf(fid,'[');
    for j = 1:1:ResourceNu
        if j < ResourceNu
            fprintf(fid,'%d,',LL(i,j));
        else
            fprintf(fid,'%d',LL(i,j));
        end
    end
    if i < length(LL)
        fprintf(fid,'],');
    else
        fprintf(fid,']');
        
        
    end
end
fprintf(fid,']');


fprintf(fid,'\n[');
for i= 1:1:ResourceNu
    if i< ResourceNu
        fprintf(fid,'%d,',max(ActPlaUseResNu(:,i)));
    else
        fprintf(fid,'%d',max(ActPlaUseResNu(:,i)));
    end
end
fprintf(fid,']');

%输出
fprintf(fid,'\n');
for i= 1:1:size(PC,1)
    fprintf(fid,'[');
    for j = 1:1:ResourceNu
        fprintf(fid,'[');
        for k=1:1:length(BB(j,:,i))
            if k<length(BB(j,:,i))
                fprintf(fid,'%d,',BB(j,k,i));
             else
                fprintf(fid,'%d',BB(j,k,i));
            end
        end
        if j < ResourceNu
            fprintf(fid,'],');
        else
            fprintf(fid,']');
        end
    end
    fprintf(fid,']\n');
end

fclose(fid);












