x = -4:0.001:12;

j = 1;
f1 = [];
for i = 1:20
    f1(i,:) = normpdf(x, j, 1);
    j = j + 0.05;
end
for i = 1:20
    y = f1';
    x = x';
    plot(x, y(:,i), '-.', 'Color', [0.65 0.85 0.43]);
    %pause(0.5);
    hold on;
end

k = 4;
f2 = [];
for i = 1:20
    f2(i,:) = normpdf(x, k, 1);
    k = k + 0.05;
end

for i = 1:20
    y2 = f2';
    x = x';
    plot(x, y2(:,i), '-.', 'Color', [1 0.4 0.6]');
    %pause(0.5);
    hold on;
end

h = 7;
f3 = [];
for i = 1:10
    f3(i,:) = normpdf(x, h, 1);
    h = h + 0.05;
end

for i = 1:10
    y3 = f3';
    x = x';
    plot(x, y3(:,i), '-.', 'Color', [0.51 0.87 1]);
    %pause(0.5);
    hold on;
end

data = [f1; f2; f3];
save('data120.mat','data');


N = size(data,1);
% epsilon = 0.01;


%
zt=data;
tic
W=[];
for i=1:size(zt,1)
    for j=1:size(zt,1)
        if i==j
            W(i,j)=0;
        else
            W(i,j)=sum((zt(i,:)-zt(j,:)).^2);
        end
    end
end

%Dieu chinh lai sai so
for i=1:length(W)
    for j=1:length(W)
        if W(i,j)<0
            W(i,j)=0;
        end
    end
end

W;
%Tinh Ws
Ws=sum(sum(W))/(length(W)*length(W)-length(W));

%Tinh ma tran K lamda
k=[];
for i=1:length(W)
    for j=1:length(W)
        if W(i,j)<Ws
            k(i,j)=exp(-W(i,j)/(Ws/5));
        else
            k(i,j)=0;
        end
    end
end

k;
%Tinh z(t+1)
ztmoi=[];
for i=1:length(W)
    tu=[];
        for l=1:length(x)
            tu=[tu,0];
        end
        mau=0;
        for j=1:length(W)
            tu=tu+zt(j,:)*k(i,j);
            mau=mau+k(i,j);
        end
    tu/mau;
    ztmoi=[ztmoi;tu/mau];
end

%So sanh chuan va do lech toi da cho phep
exilanh=10^-2;
max(max(abs(ztmoi-zt)));
vonglap=0;
while max(max(abs(ztmoi-zt)))>exilanh 
    vonglap=vonglap+1
    zt=ztmoi;
    figure
    
% for i=1:N
%      da(i,:)=plot(x,data(i,:),'-.b');
%      set(gcf,'color','w');
%      hold on
%      mn(i,:)=plot(x,zt(i,:),'-r');
%      set(gcf,'color','w');
% end

% for i=1:N
%      da(i,:)=plot(x,data(i,:),'-.b');
%      set(gcf,'color','w');
%      hold on
% end
% 
% for i=1:N
% mn(i,:)=plot(x,zt(i,:),'-r');
%      set(gcf,'color','w');
%      hold on
% end


for i = 1:N
    da(i,:) = plot(x, data(i,:), '-.', 'Color', [0.2 0.6 1]); % Màu tím nh?t
    set(gcf, 'color', 'w');
    hold on
end


for i = 1:N
    mn(i,:) = plot(x, zt(i,:), '-.', 'Color', [1 0.23 0.47]); % Màu ?? nh?t
    set(gcf, 'color', 'w');
    hold on
end



legend([da(1,:) mn(1,:)],{'Input PDF','Output PDF'});
 hold off
 
 W=[];
for i=1:size(zt,1)
    for j=1:size(zt,1)
        if i==j
            W(i,j)=0;
        else
        W(i,j)=sum((zt(i,:)-zt(j,:)).^2);
        end
    end
end

%Dieu chinh lai sai so
for i=1:length(W)
    for j=1:length(W)
        if W(i,j)<0
            W(i,j)=0;
        end
    end
end

W;
% alp=alpha1(alp,k);
%Tinh ma tran K lamda
k=[];
for i=1:length(W)
    for j=1:length(W)
        if W(i,j)<Ws
            k(i,j)=exp(-W(i,j)/(Ws/5));
        else
            k(i,j)=0;
        end
    end
end

k;
%Tinh z(t+1)
ztmoi=[];
for i=1:length(W)
    tu=[];
        for l=1:length(x)
        tu=[tu,0];
        end
        mau=0;
        for j=1:length(W)
        tu=tu+zt(j,:)*k(i,j);
        mau=mau+k(i,j);
        end
    tu/mau;
    ztmoi=[ztmoi;tu/mau];
end
max(max(abs(ztmoi-zt)));
end
%

f_new = ztmoi;


figure

% for i=1:N
%     da(i,:) = plot(x,data(i,:),'-.b');
%     set(gcf,'color','w');
%     hold on
%     %pause(0.5);
% end
% 
% for i=1:N
%     mn(i,:) = plot(x,f_new(i,:),'-.r');
%     set(gcf,'color','w');
%     %pause(0.5);
% end


for i=1:N
    da(i,:) = plot(x,data(i,:),'-.', 'Color', [0.2 0.6 1]);
    set(gcf,'color','w');
    hold on
    %pause(0.5);
end

for i=1:N
    mn(i,:) = plot(x,f_new(i,:),'-.', 'Color', [1 0.23 0.47]);
    set(gcf,'color','w');
    %pause(0.5);
end


legend([da(1,:) mn(1,:)],{'Input PDF','Output PDF'});
hold off

Khoang_cach = [];
for j=1:N
    for i=1:N
        Khoang_cach(i,j) = sum((zt(i,:)-zt(j,:)).^2);;
    end
end

% Chon u phu hop de phan tach cac chum
% u = 1;
u = min(Khoang_cach(find(Khoang_cach>0.2)))-0.1;


for i=1:N
    for j=1:N
        if Khoang_cach(i,j) < u
           E(i,j) = j;
        else E(i,j) = 0;
        end
    end
end

H = unique(E, 'rows');

Q = flipud(H);

for i=1:size(Q,1)
    L(find(Q(i,:)~=0)) = i;
end

unique_values = unique(L);
% So chum
num_unique = numel(unique_values);
so_chum = num_unique;

disp(L);
disp(so_chum);




U = zeros(so_chum,N);

for j=1:N
    for i=1:so_chum
        if L(j) == i
            U(i,j) = 1;
        else U(i,j) =0;
        end
    end
end

f_bar_m = sum(U,2);
f_bar_t = zeros(size(U,1),size(data,2));

for i=1:size(U,1)
    for j=1:size(U,2)
        ff_bar = U(i,j)*data(j,:);
        f_bar_t(i,:) = f_bar_t(i,:) + ff_bar;
    end
end

f_bar = zeros(size(f_bar_t));
for i=1:size(f_bar_t,1)
    f_bar(i,:) = f_bar_t(i,:)/f_bar_m(i,:);
end

U_new = U;
num_steps = 0;

epsilon = 0.0001;

while true
    U_old = U_new;
    for j = 1:size(U,2)
        for i = 1:size(U,1)
            sum_val = 0; 
            if norm(data(j,:)-f_bar(i,:), 2) == 0
                U_new(i, j) = 1;
            else
                for m = 1:size(U,1)
                    sum_val = sum_val + (sum((data(j,:)-f_bar(i,:)).^2) / ...
                        sum((data(j,:)-f_bar(m,:))).^2);
                end
            U_new(i, j) = 1 / sum_val;
            end
        end
        %U_new(setdiff(1:size(U_new, 1), find(U_new(:, j) == 1)), j) = 0;
    end
    num_steps = num_steps + 1; % Tong so buoc lap sau moi lan lap
    
            f_bar_m = sum(U_new,2);

            f_bar_t = zeros(size(U_new,1),size(data,2));
            
            for i=1:size(U,1)
                for j=1:size(U,2)
                    ff_bar = U(i,j)*data(j,:);
                    f_bar_t(i,:) = f_bar_t(i,:) + ff_bar;
                end
            end
            
            f_bar = zeros(size(f_bar_t));
            for i=1:size(f_bar_t,1)
                f_bar(i,:) = f_bar_t(i,:)/f_bar_m(i,:);
            end
    
    if max(abs(U_new - U_old)) < epsilon
        break;
    end
end

%disp(num_steps);
disp(U_new);


