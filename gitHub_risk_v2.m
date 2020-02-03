%% Tthe probability of introducing at least one case from Wuhan to city j by day t
%  Written by Dr. Zhanwei Du from The University of Texas at Austin, Austin, Texas 78712, The United States of America
%  Email: duzhanwei0@gmail.com

%% Input
% Pop_wuhan: Population size of Wuhan city
% Iwt:the number of infected cases at time t
% cityNames: list of city names
% wuhanID: the index of Wuhan in cityNames
% Jan23: Timing of Jan.23, 2020
% Tpeirod: The length of studied period
% Mobs: The mobility matrix between cities
%% Output
% prob_ToFrom: the probability of introducing at least one case from Wuhan to city j by day t

%% To
pi = Iwt/Pop_wuhan;
prob_To = zeros(Tpeirod,size(Mobs,1));
for i=1:length(cityNames)
    for t=1:Tpeirod
        temp = 0;
        for j=1:t
            temp = temp+Mobs(i,wuhanID,j) *pi(j);
        end
        prob_To(t,i) = 1-exp(-temp) ;
    end
end

%% From
prob_From = zeros(Tpeirod,size(Mobs,1));
for i=1:length(cityNames)
    for t=1:Tpeirod
        temp = 0;
        for j=1:t
            temp = temp+Mobs(wuhanID,i,j) *pi(j);
        end
        prob_From(t,i) = 1-exp(-temp) ;
    end
end

%%
prob_ToFrom = zeros( Tpeirod,length(cityNames));
for i=1:length(cityNames)
    for t=1:Tpeirod
        prob_ToFrom(t,i) = 1-(1-prob_To(t,i))*(1-prob_From(t,i));
    end
end


